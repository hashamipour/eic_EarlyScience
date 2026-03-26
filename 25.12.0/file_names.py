import os
import sys
from pathlib import Path

folder = sys.argv[1] if len(sys.argv) > 1 else "./"

folder_path = Path(folder).resolve()
folder_name = folder_path.name
output_file = f"{folder_name}.txt"

root_files = [str(folder_path / f) for f in os.listdir(folder) if f.endswith(".root")]

with open(output_file, "w") as f:
    for rf in root_files:
        f.write(rf + "\n")

print(f"Found {len(root_files)} .root files, saved to {output_file}")