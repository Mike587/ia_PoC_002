"""
czi_analysis.py — ANALYSIS-SPECIFIC logic for ia_PoC_002 (nucleus detection).

This is the part you REPLACE when building a new SmartMic image analysis. None
of it is part of the SmartMic contract (that lives in `analyze_czi.py`). Here
that contract is satisfied by:

  * get_scene_center_positions() — read the per-scene stage centre from the CZI
    metadata, so detected features can be expressed as absolute stage coords.
  * find_nuclei()                — the actual detection; returns objects carrying
    `abs_x_m` / `abs_y_m` (the only fields SmartMic strictly needs).
  * save_result_image()          — optional QC overlay (not consumed by SmartMic).

Swap these for your own detection while keeping the absolute-coordinate output.
"""
from pathlib import Path

import numpy as np
from scipy import ndimage as ndi
from skimage.filters import threshold_otsu, gaussian
from skimage.morphology import opening, disk
from skimage.segmentation import watershed, clear_border
from skimage.feature import peak_local_max
from skimage.measure import regionprops, label
from bioio import BioImage


def _read_actual_axis_positions(root) -> tuple[float | None, float | None]:
    """
    Read actual stage encoder positions (µm) recorded at acquisition time
    from MTBStageAxisX / MTBStageAxisY ParameterCollection entries,
    and convert to meters. Falls back to None if not found.
    """
    x_m = y_m = None
    for pc in root.iter("ParameterCollection"):
        cid = pc.get("Id", "")
        pos_elem = pc.find("Position")
        if pos_elem is None or pos_elem.text is None:
            continue
        try:
            val_m = float(pos_elem.text) / 1e6
        except ValueError:
            continue
        if cid == "MTBStageAxisX":
            x_m = val_m
        elif cid == "MTBStageAxisY":
            y_m = val_m
    return x_m, y_m


def get_scene_center_positions(img: BioImage) -> list[dict]:
    """
    Return stage center positions (m) for each scene.
    Prefers actual hardware encoder positions (MTBAxisX/Y) over the
    planned CenterPosition, as the stage may not land exactly on target.
    """
    root = img.metadata
    actual_x_m, actual_y_m = _read_actual_axis_positions(root)

    positions = []
    for scene in root.iter("Scene"):
        cp = scene.findtext("CenterPosition")
        if cp:
            planned_x_m, planned_y_m = (float(v) / 1e6 for v in cp.split(","))
        else:
            planned_x_m = planned_y_m = None

        # Use planned CenterPosition as the stage reference.
        # The encoder readback (MTBStageAxisX/Y) has IsPrecise=false and
        # overshoots more than the nominal position, so planned is more accurate.
        center_x_m = planned_x_m
        center_y_m = planned_y_m

        positions.append({
            "index": scene.get("Index"),
            "name": scene.get("Name"),
            "center_x_m": center_x_m,
            "center_y_m": center_y_m,
            "planned_x_m": planned_x_m,
            "planned_y_m": planned_y_m,
            "actual_x_m": actual_x_m,
            "actual_y_m": actual_y_m,
        })
    return positions


def find_nuclei(img: BioImage, scene_center: dict, max_nuclei: int = 1000,
                min_area_m2: float = 25e-12):
    """
    Detect nuclei in the DAPI channel via Otsu threshold + watershed.
    Rejects nuclei touching the image border and smaller than min_area_m2.
    Returns (nuclei list, raw array, label image, centroid array [row, col]).
    """
    pixel_size_y = img.physical_pixel_sizes.Y * 1e-6  # µm -> m
    pixel_size_x = img.physical_pixel_sizes.X * 1e-6

    arr = img.get_image_data("YX", T=0, C=0, Z=0).astype(np.float32)
    height, width = arr.shape

    smoothed = gaussian(arr, sigma=2)
    thresh = threshold_otsu(smoothed)
    binary = smoothed > thresh
    binary = opening(binary, disk(3))

    distance = ndi.distance_transform_edt(binary)
    peak_coords = peak_local_max(distance, min_distance=15, labels=binary)
    peak_mask = np.zeros_like(distance, dtype=bool)
    peak_mask[tuple(peak_coords.T)] = True
    markers = label(peak_mask)
    segmented = watershed(-distance, markers, mask=binary)
    segmented = clear_border(segmented)

    props = regionprops(segmented)
    cx_m = scene_center["center_x_m"]
    cy_m = scene_center["center_y_m"]

    nuclei = []
    centroids_rc = []
    for region in props:
        area_m2 = region.area * pixel_size_x * pixel_size_y
        if area_m2 < min_area_m2:
            continue
        if len(nuclei) >= max_nuclei:
            break
        row, col = region.centroid
        dx = (col - width / 2) * pixel_size_x
        dy = (row - height / 2) * pixel_size_y
        nuclei.append({
            "id": int(region.label),
            "abs_x_m": cx_m + dx,
            "abs_y_m": cy_m + dy,
            "area_m2": area_m2,
            "centroid_col": int(round(col)),
            "centroid_row": int(round(row)),
        })
        centroids_rc.append([row, col])

    return nuclei, arr, segmented, np.array(centroids_rc)


def save_result_image(arr, nuclei, centroids_rc, output_path: Path):
    """Save a DAPI image with nucleus IDs annotated at each centroid."""
    import matplotlib
    matplotlib.use("Agg")  # non-interactive backend, safe for scripts
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(16, 12), dpi=150)

    # Display with percentile stretch so dim nuclei are visible
    vmin, vmax = np.percentile(arr, (1, 99.5))
    ax.imshow(arr, cmap="gray", vmin=vmin, vmax=vmax, interpolation="nearest")

    for nucleus, (row, col) in zip(nuclei, centroids_rc):
        ax.text(
            col, row,
            str(nucleus["id"]),
            fontsize=8,
            color="orange",
            ha="center",
            va="center",
            fontweight="bold",
        )

    ax.set_title(f"Nucleus detection — {len(nuclei)} nuclei (edge-touching excluded)")
    ax.axis("off")
    plt.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)
