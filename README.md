# ia_PoC_002 — image-analysis scaffold for SmartMic

A standalone CZI image-analysis project that runs in its **own pixi environment**,
so its analysis dependencies (bioio, scikit-image, scipy, matplotlib) stay
isolated from the ZEN-API environment used to drive the microscope.

**This is an example / scaffold, not a fixed tool.** It is meant as a template
for *any* kind of image analysis you want to plug into the SmartMic pipeline
([Mike587/SmartMic](https://github.com/Mike587/SmartMic)). This particular PoC
happens to do DAPI nucleus detection, but the point is the integration pattern:
SmartMic hands an analysis script a CZI and reads back a JSON of target
positions. To build a different analysis (organoids, colonies, debris rejection,
custom features, …), copy this project and swap the detection logic while
keeping the same input/output contract described below.

## How SmartMic uses it

SmartMic invokes this project as a **subprocess** from
`MS_SmartMic_PoC.py` (via `MS_image_analysis.run_analysis`), running it in *this*
project's pixi environment. The contract is deliberately simple and
language/-dependency-agnostic:

- **Input:** SmartMic passes a CZI path + an output directory (and optional
  `--prefix` / `--offset`) on the command line.
- **Output:** the script writes a `*_nuclei.json` into the output directory;
  SmartMic reads that JSON back and uses the absolute stage coordinates to drive
  the stage to each detected target for detailed acquisition.

Keep that contract and SmartMic does not care what happens in between — that is
what makes this a reusable scaffold.

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

## Output

Each run writes three files to `<output_dir>`, all sharing the optional
`{prefix}_` stem:

| File | Description |
|------|-------------|
| `{prefix}_nuclei.json` | **The machine-readable result SmartMic consumes** — a JSON array of detected targets (see spec below). |
| `{prefix}_analysis.log` | Human-readable log of the run (image shape, pixel size, per-scene stage positions, detection summary). |
| `{prefix}_nuclei_overlay.png` | DAPI image with each detected nucleus ID annotated at its centroid — for visual QC. |

### `*_nuclei.json` specification

A JSON **array** of objects, one per detected target. Edge-touching regions and
regions smaller than `min_area_m2` (default 25 µm²) are excluded; at most
`max_nuclei` (default 1000) are returned. An empty array (`[]`) means nothing
was detected.

| Field | Type | Units | Meaning |
|-------|------|-------|---------|
| `id` | int | — | Region label, unique within this image. |
| `abs_x_m` | float | metres | **Absolute stage X** of the target centroid (scene centre + in-image offset, plus any `--offset`). This is what SmartMic moves the stage to. |
| `abs_y_m` | float | metres | **Absolute stage Y** of the target centroid. |
| `area_m2` | float | m² | Segmented area of the region. |
| `centroid_col` | int | pixels | Centroid column (X) in the image — used for the overlay annotation. |
| `centroid_row` | int | pixels | Centroid row (Y) in the image. |

> Note on the contract: SmartMic only requires the absolute stage coordinates
> (`abs_x_m`, `abs_y_m`) to revisit each target. The other fields are useful
> metadata/QC. A different analysis can add its own fields freely, as long as it
> keeps `abs_x_m` / `abs_y_m` so the pipeline can navigate to each hit.

### Example

```json
[
  {
    "id": 1,
    "abs_x_m": 0.0441234,
    "abs_y_m": 0.0093120,
    "area_m2": 6.8421e-11,
    "centroid_col": 512,
    "centroid_row": 387
  },
  {
    "id": 4,
    "abs_x_m": 0.0442871,
    "abs_y_m": 0.0091044,
    "area_m2": 4.1130e-11,
    "centroid_col": 901,
    "centroid_row": 145
  }
]
```
