import argparse
import xml.etree.ElementTree as ET
from xml.dom import minidom
from bioio import BioImage
from czi_analysis import get_scene_center_positions, find_nuclei

parser = argparse.ArgumentParser(description="Explore CZI metadata and detect nuclei.")
parser.add_argument("czi_file", nargs="?", default="experiment_pos_000077.czi")
args = parser.parse_args()

img = BioImage(args.czi_file)

print(f"Shape: {img.shape}")
print(f"Dims:  {img.dims}")
print(f"Physical pixel sizes: {img.physical_pixel_sizes}")
print(f"Channel names: {img.channel_names}")
print(f"Scene count: {len(img.scenes)}")
print(f"Current scene: {img.current_scene}")

print("\n--- Raw metadata (XML) ---")
meta_xml = ET.tostring(img.metadata, encoding="unicode")
print(minidom.parseString(meta_xml).toprettyxml(indent="  "))

positions = get_scene_center_positions(img)

print("\n--- Scene center positions ---")
for p in positions:
    print(f"  Scene {p['index']} ({p['name']})")
    print(f"    Actual (encoder): X={p['actual_x_m']} m, Y={p['actual_y_m']} m")
    print(f"    Planned:          X={p['planned_x_m']} m, Y={p['planned_y_m']} m")

print("\n--- Nucleus detection ---")
nuclei, *_ = find_nuclei(img, scene_center=positions[0])
print(f"Found {len(nuclei)} nuclei (max 1000, edge-touching excluded)\n")
print(f"{'ID':>4}  {'Abs X (m)':>18}  {'Abs Y (m)':>18}  {'Area (m²)':>14}")
print("-" * 62)
for n in nuclei:
    print(f"{n['id']:>4}  {n['abs_x_m']:>18.6e}  {n['abs_y_m']:>18.6e}  {n['area_m2']:>14.4e}")
