#!/usr/bin/env python3
"""
Script 3: Calculate Distances (Zero-RAM Streaming Version)
"""
import subprocess
import sys
import argparse
import time
from pathlib import Path


def main():
    # 1. Argument Parsing
    parser = argparse.ArgumentParser(description="Step 3: Calculate Distances")
    parser.add_argument("input_dir", help="Directory with .msh file")
    parser.add_argument("-o", "--output-dir", required=True)
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()

    # 2. Setup Paths
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    input_dir = Path(args.input_dir)

    sketch_path = input_dir / "genomes_sketch.msh"
    dist_file = out_dir / "distances.txt"
    dist_file_tmp = out_dir / "distances.txt.tmp"

    if not sketch_path.exists():
        print(f"[ERROR] Sketch file not found: {sketch_path}", file=sys.stderr)
        sys.exit(1)

    # 3. Define Mash Command
    cmd = [
        "mash", "dist",
        "-t",                    # Tabular output (square matrix)
        "-p", str(args.threads),
        str(sketch_path),        # Reference set
        str(sketch_path)         # Query set — triggers All-vs-All
    ]

    print("-" * 60)
    print(f"[INFO] Starting calculation (All-vs-All)")
    print(f"[INFO] Threads:      {args.threads}")
    print(f"[INFO] Output:       {dist_file}")
    print(f"[INFO] Temp output:  {dist_file_tmp}")
    print("-" * 60)

    # Clean up any leftover tmp file from a previous failed run
    if dist_file_tmp.exists():
        print(f"[WARNING] Stale temp file found, removing: {dist_file_tmp}")
        dist_file_tmp.unlink()

    # 4. Execution via Direct Streaming
    start_time = time.time()
    try:
        with open(dist_file_tmp, "w") as f_out:
            subprocess.run(
                cmd,
                stdout=f_out,
                stderr=subprocess.PIPE,
                text=True,
                check=True
            )

        # 5. Atomic rename: only happens if subprocess did not raise
        dist_file_tmp.rename(dist_file)

        duration = time.time() - start_time
        print(f"[SUCCESS] Calculation finished in {duration / 60:.2f} minutes.")
        print(f"[SUCCESS] Matrix saved to {dist_file}")

    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Mash dist failed (exit code {e.returncode}).", file=sys.stderr)
        print(f"[ERROR] Stderr: {e.stderr}", file=sys.stderr)
        if dist_file_tmp.exists():
            dist_file_tmp.unlink()
            print(f"[INFO] Partial temp file removed.", file=sys.stderr)
        sys.exit(1)

    except Exception as e:
        print(f"[ERROR] Unexpected error: {e}", file=sys.stderr)
        if dist_file_tmp.exists():
            dist_file_tmp.unlink()
            print(f"[INFO] Partial temp file removed.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()