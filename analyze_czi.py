"""
analyze_czi.py — entry point for a SmartMic image-analysis project.

This file is the SCAFFOLD shared by every analysis plugged into the SmartMic
pipeline.  SmartMic (MS_image_analysis.run_analysis) always launches an analysis
as a subprocess with this fixed CLI:

    pixi run python analyze_czi.py <path_to_czi> <output_dir> --prefix <tag>

and then reads back the `<prefix>_targets.json` written into <output_dir>.

The bits marked "GENERAL" below are the contract with SmartMic and should be
kept (largely unchanged) in any new analysis.  The bit marked
"ANALYSIS-SPECIFIC" is what you replace to do a different analysis — see
`czi_analysis.py`.

Writes to <output_dir>/:
  {prefix}_analysis.log        — human-readable log of the run        (GENERAL)
  {prefix}_targets.json        — detected targets as a JSON array     (result SmartMic reads)
  {prefix}_nuclei_overlay.png  — QC overlay image                     (analysis-specific)

If --prefix is omitted the filenames have no prefix.
"""
import argparse
import json
import logging
import sys
from datetime import datetime
from pathlib import Path

from bioio import BioImage
from czi_analysis import get_scene_center_positions, find_nuclei, save_result_image


# ============================ GENERAL ============================
# Logging setup: file + stdout.  SmartMic captures stdout into the run log,
# so logging to stdout is part of how results/diagnostics flow back.  Keep this
# in any analysis.
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
    # ======================== GENERAL ========================
    # The SmartMic analysis contract: positional <image> <output> plus
    # --prefix <tag>.  Do NOT change this CLI signature — MS_image_analysis
    # depends on it so any analysis project is interchangeable.
    parser = argparse.ArgumentParser(description="Detect nuclei in a CZI DAPI image.")
    parser.add_argument("image", help="Path to the .czi file")
    parser.add_argument("output", help="Output directory (created if absent)")
    parser.add_argument(
        "--prefix", default="",
        help="Filename prefix for all output files (e.g. D3_P1 → D3_P1_targets.json)",
    )
    args = parser.parse_args()

    image_path = Path(args.image)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Output paths share the {prefix}_ stem.  `json_path` is the result file
    # SmartMic reads back — every analysis must write its targets here.
    stem = f"{args.prefix}_" if args.prefix else ""
    log_path  = output_dir / f"{stem}analysis.log"
    json_path = output_dir / f"{stem}targets.json"

    log = setup_logging(log_path)
    log.info(f"Started analysis  file={image_path}  output={output_dir}")
    log.info(f"Timestamp: {datetime.now().isoformat()}")

    # Fail fast with a clear log + non-zero exit if the input is missing.
    # SmartMic treats a non-zero exit (or a missing result JSON) as "no targets".
    if not image_path.exists():
        log.error(f"Image not found: {image_path}")
        sys.exit(1)

    # ==================== ANALYSIS-SPECIFIC ====================
    # Everything from here until the result-write is what you replace for a
    # different analysis.  The ONLY hard requirement is: write the detected
    # targets to `json_path` as a JSON array of objects that include absolute
    # stage coordinates `abs_x_m` / `abs_y_m` (the fields SmartMic uses to drive
    # the stage), and exit 0 on success.
    log.info("Loading image...")
    img = BioImage(str(image_path))
    log.info(f"Shape:              {img.shape}")
    log.info(f"Dims:               {img.dims}")
    log.info(f"Channel names:      {img.channel_names}")
    log.info(f"Physical pixel size Y={img.physical_pixel_sizes.Y} µm  X={img.physical_pixel_sizes.X} µm")
    log.info(f"Scene count:        {len(img.scenes)}")

    positions = get_scene_center_positions(img)
    if not positions:
        log.error("No <Scene> elements found in CZI metadata; cannot determine stage positions.")
        sys.exit(1)
    for p in positions:
        def fmt(v): return f"{v:.6e} m" if v is not None else "n/a"
        log.info(f"Scene {p['index']} ({p['name']}): "
                 f"actual X={fmt(p['actual_x_m'])}  Y={fmt(p['actual_y_m'])}  "
                 f"(planned X={fmt(p['planned_x_m'])}  Y={fmt(p['planned_y_m'])})")

    max_nuclei = 1000
    log.info("Running nucleus detection (Otsu + watershed, edge-touching rejected, min area 25 µm²)...")
    nuclei, _arr, _seg, _cents, truncated = find_nuclei(img, scene_center=positions[0], max_nuclei=max_nuclei)
    log.info(f"Detected {len(nuclei)} nuclei (edge-touching excluded, max {max_nuclei})")
    if truncated:
        log.warning(f"Nucleus count reached max_nuclei={max_nuclei}; "
                     "some qualifying nuclei were left out of the output.")

    if nuclei:
        areas = [n["area_m2"] for n in nuclei]
        log.info(f"Area min={min(areas):.4e} m²  max={max(areas):.4e} m²  mean={sum(areas)/len(areas):.4e} m²")

    # ======================== GENERAL ========================
    # Write the result JSON SmartMic reads back.  The write mechanism is
    # general; the per-object schema is analysis-specific (see README).
    log.info(f"Writing results to {json_path}")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(nuclei, f, indent=2)

    # ==================== ANALYSIS-SPECIFIC ====================
    # Optional QC overlay — not consumed by SmartMic, purely for humans.
    img_path = output_dir / f"{stem}nuclei_overlay.png"
    log.info(f"Saving result image to {img_path}")
    save_result_image(_arr, nuclei, _cents, img_path)

    log.info("Done.")


if __name__ == "__main__":
    main()
