"""
Usage: python analyze_czi.py <path_to_image> <path_to_output> [--offset DX DY] [--prefix PREFIX]

Writes to <path_to_output>/:
  {prefix}_analysis.log      — human-readable log of the analysis run
  {prefix}_nuclei.json       — detected nuclei as a JSON array
  {prefix}_nuclei_overlay.png

If --prefix is omitted the filenames have no prefix (original behaviour).

Optional offset correction (in metres) for the camera centre / optical axis
calibration not stored in the CZI.  Example — 44 µm in X, 9 µm in Y:
  --offset 44e-6 9e-6
"""
import argparse
import json
import logging
import sys
from datetime import datetime
from pathlib import Path

from bioio import BioImage
from czi_analysis import get_scene_center_positions, find_nuclei, save_result_image


def setup_logging(log_path: Path) -> logging.Logger:
    logger = logging.getLogger("analyze_czi")
    logger.setLevel(logging.DEBUG)
    fmt = logging.Formatter("%(asctime)s  %(levelname)-8s  %(message)s", datefmt="%Y-%m-%dT%H:%M:%S")
    logger.addHandler(logging.FileHandler(log_path, encoding="utf-8"))
    logger.addHandler(logging.StreamHandler(sys.stdout))
    for h in logger.handlers:
        h.setFormatter(fmt)
    return logger


def main():
    parser = argparse.ArgumentParser(description="Detect nuclei in a CZI DAPI image.")
    parser.add_argument("image", help="Path to the .czi file")
    parser.add_argument("output", help="Output directory (created if absent)")
    parser.add_argument(
        "--offset", nargs=2, type=float, default=[0.0, 0.0],
        metavar=("DX", "DY"),
        help="Camera centre offset correction in metres (e.g. 44e-6 9e-6)",
    )
    parser.add_argument(
        "--prefix", default="",
        help="Filename prefix for all output files (e.g. D3_P1 → D3_P1_nuclei.json)",
    )
    args = parser.parse_args()
    offset_x_m, offset_y_m = args.offset

    image_path = Path(args.image)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    stem = f"{args.prefix}_" if args.prefix else ""
    log_path  = output_dir / f"{stem}analysis.log"
    json_path = output_dir / f"{stem}nuclei.json"

    log = setup_logging(log_path)
    log.info(f"Started analysis  file={image_path}  output={output_dir}")
    log.info(f"Timestamp: {datetime.now().isoformat()}")

    if not image_path.exists():
        log.error(f"Image not found: {image_path}")
        sys.exit(1)

    log.info("Loading image...")
    img = BioImage(str(image_path))
    log.info(f"Shape:              {img.shape}")
    log.info(f"Dims:               {img.dims}")
    log.info(f"Channel names:      {img.channel_names}")
    log.info(f"Physical pixel size Y={img.physical_pixel_sizes.Y} µm  X={img.physical_pixel_sizes.X} µm")
    log.info(f"Scene count:        {len(img.scenes)}")

    positions = get_scene_center_positions(img)
    for p in positions:
        def fmt(v): return f"{v:.6e} m" if v is not None else "n/a"
        log.info(f"Scene {p['index']} ({p['name']}): "
                 f"actual X={fmt(p['actual_x_m'])}  Y={fmt(p['actual_y_m'])}  "
                 f"(planned X={fmt(p['planned_x_m'])}  Y={fmt(p['planned_y_m'])})")

    if offset_x_m or offset_y_m:
        log.info(f"Applying camera centre offset: DX={offset_x_m:.4e} m  DY={offset_y_m:.4e} m")
        positions[0]["center_x_m"] += offset_x_m
        positions[0]["center_y_m"] += offset_y_m

    log.info("Running nucleus detection (Otsu + watershed, edge-touching rejected, min area 25 µm²)...")
    nuclei, _arr, _seg, _cents = find_nuclei(img, scene_center=positions[0], max_nuclei=1000)
    log.info(f"Detected {len(nuclei)} nuclei (edge-touching excluded, max 1000)")

    if nuclei:
        areas = [n["area_m2"] for n in nuclei]
        log.info(f"Area min={min(areas):.4e} m²  max={max(areas):.4e} m²  mean={sum(areas)/len(areas):.4e} m²")

    log.info(f"Writing results to {json_path}")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(nuclei, f, indent=2)

    img_path = output_dir / f"{stem}nuclei_overlay.png"
    log.info(f"Saving result image to {img_path}")
    save_result_image(_arr, nuclei, _cents, img_path)

    log.info("Done.")


if __name__ == "__main__":
    main()
