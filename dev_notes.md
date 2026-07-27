# ia_PoC_002 — dev notes (learnings + decisions)

Empirical findings that aren't visible from reading a single function — mirrors
the format of SmartMic's own `DEV_NOTES.md`. Kept in this repo because
`ia_PoC_002` is a standalone project with its own git history.

---

## Pixel → stage axis mapping (`czi_analysis.find_nuclei`)

`find_nuclei` converts a detected nucleus centroid to absolute stage
coordinates as `abs_x_m = center_x_m + (col - width/2) * pixel_size_x` and the
Y equivalent for rows. This assumes image column/row map directly onto stage
X/Y with no flip or transpose.

**Confirmed on hardware (per project owner, 2026-07-27): this sign/orientation
convention is correct** — `abs_x_m`/`abs_y_m` drives SmartMic back to the same
physical feature it was computed from. Not independently re-verified as part
of this note; recorded here so the assumption isn't silently re-litigated or
mistaken for unverified in a future review. If the camera orientation,
binning, or objective/optovar path ever changes on this scope, re-check this
before trusting `abs_x_m`/`abs_y_m` again.

## Scene center: planned position preferred over encoder readback

`get_scene_center_positions` uses the CZI metadata's planned `<CenterPosition>`
as the stage reference, not the `MTBStageAxisX`/`MTBStageAxisY` encoder
readback (also read and returned, but only for logging/QC). Reason: the
encoder readback has `IsPrecise=false` and overshoots more than the nominal
(planned) position, so planned is the more accurate reference for computing
absolute nucleus coordinates. See the comment at `get_scene_center_positions`
in `czi_analysis.py`.
