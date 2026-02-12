#!/usr/bin/env python3
"""
Script 2: Calculate Distances
"""
import subprocess
import sys
import argparse
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description="Step 2: Calculate Distances")
    parser.add_argument("input_dir", help="Directory with .msh file")
    parser.add_argument("-o", "--output-dir", required=True)
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    sketch_path = Path(args.input_dir) / "genomes_sketch.msh"
    if not sketch_path.exists():
        print(f"[ERROR] Sketch file not found: {sketch_path}")
        sys.exit(1)

    dist_file = out_dir / "distances.txt"

    cmd = [
        "mash", "dist",
        "-t", # tabular output
        "-p", str(args.threads),
        str(sketch_path),
        str(sketch_path)
    ]

    print("[INFO] Calculating distances (All vs All)...")
    try:
        with open(dist_file, "w") as f:
            subprocess.run(cmd, stdout=f, check=True)
        print(f"[SUCCESS] Distances saved to {dist_file}")
    except subprocess.CalledProcessError:
        print("[ERROR] Mash dist failed.")
        sys.exit(1)

if __name__ == "__main__":
    main()