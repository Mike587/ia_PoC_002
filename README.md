# ia_PoC_002 — nuclei detection for SmartMic

Standalone CZI nucleus-detection analysis, run in its **own pixi environment** so
its image-analysis dependencies (bioio, scikit-image, scipy, matplotlib) stay
isolated from the ZEN-API environment.

This project is invoked as a subprocess by the SmartMic PoC pipeline
(`MS_SmartMic_PoC.py`, via `MS_image_analysis.run_analysis`), which passes it an
overview CZI and reads back the `*_nuclei.json` it writes.

## Files

| File | Role |
|------|------|
| `analyze_czi.py` | Entry point. Detects nuclei in a CZI and writes `*_nuclei.json`, `*_analysis.log`, `*_nuclei_overlay.png`. |
| `czi_analysis.py` | Detection logic: scene-center positions from metadata + Otsu/watershed nucleus segmentation + overlay rendering. |
| `read_czi_metadata.py` | Dev helper: dump a CZI's metadata XML and print the detected-nuclei table. |

## Usage

```sh
pixi run python analyze_czi.py <path_to_czi> <output_dir> [--prefix PREFIX] [--offset DX DY]
```

Example:

```sh
pixi run python analyze_czi.py experiment_pos_000077.czi output_test/ --prefix D3_P1
```

`--offset DX DY` applies a camera-centre offset correction in metres (e.g.
`--offset 44e-6 9e-6`) for calibration not stored in the CZI.

Inspect a CZI's metadata / detection table without writing outputs:

```sh
pixi run python read_czi_metadata.py <path_to_czi>
```
