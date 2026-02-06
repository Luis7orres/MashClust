#!/usr/bin/env python3
"""
Script 1: Sketch and Filter (Target vs Non-Target)
Updated to support 'Cluster All' mode via --no-filter.
"""
import subprocess
import sys
import argparse
from pathlib import Path

def main():
    # 1. Argument Parsing
    parser = argparse.ArgumentParser(description="Step 1: Sketch and Filter")
    parser.add_argument("genomes_dir", help="Directory containing genome folders from Step 0")
    parser.add_argument("-o", "--output-dir", required=True)
    # Changed --filter to not required if --no-filter is present
    parser.add_argument("-f", "--filter", help="Bacteria name prefix (e.g., staphylococcus_aureus)")
    parser.add_argument("--no-filter", action="store_true", help="Disable filtering and sketch all genomes")
    parser.add_argument("-k", "--kmer", type=int, default=31)
    parser.add_argument("-s", "--sketch-size", type=int, default=100000)
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()

    # Validation: Must have either a filter or the no-filter flag
    if not args.no_filter and not args.filter:
        print("[ERROR] You must provide a filter (-f) or use the --no-filter flag.")
        sys.exit(1)

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    target_genomes = []
    non_target_genomes = []
    
    print(f"[INFO] Scanning genomes in {args.genomes_dir}")
    if args.no_filter:
        print(f"[INFO] Filtering DISABLED: All genomes will be treated as targets.")
    else:
        print(f"[INFO] Filtering for folders starting with: '{args.filter}'")

    base_path = Path(args.genomes_dir)
    
    # 2. Categorizing Genomes
    # We search for .fna files recursively
    for genome_file in base_path.rglob("GENOME_*.fna"):
        # Parent folder name (e.g., staphylococcus_aureus_GCF...)
        parent_folder = genome_file.parent.name.lower()
        
        if args.no_filter:
            # Mode: Cluster everything
            target_genomes.append(genome_file)
        elif args.filter and parent_folder.startswith(args.filter.lower()):
            # Mode: Filter by name
            target_genomes.append(genome_file)
        else:
            # Mode: Not a target
            non_target_genomes.append(genome_file)

    print(f"[INFO] Targets found (to be sketched): {len(target_genomes)}")
    print(f"[INFO] Non-targets found (skipped):     {len(non_target_genomes)}")

    # 3. Save lists for later steps (Required for Snakemake and Script 5)
    with open(out_dir / "targets.txt", "w") as f:
        for p in target_genomes:
            f.write(f"{p}\n")
            
    with open(out_dir / "non_targets.txt", "w") as f:
        for p in non_target_genomes:
            f.write(f"{p}\n")

    if not target_genomes:
        print("[ERROR] No target genomes found! Check your input directory or filter.")
        sys.exit(1)

    # 4. Mash Sketch Execution
    sketch_file = out_dir / "genomes_sketch"
    
    cmd = [
        "mash", "sketch",
        "-s", str(args.sketch_size),
        "-k", str(args.kmer),
        "-p", str(args.threads),
        "-o", str(sketch_file),
        "-l", str(out_dir / "targets.txt") # Use the list we just created
    ]

    print("-" * 60)
    print(f"[INFO] Starting Mash Sketch on {len(target_genomes)} genomes...")
    print(f"[INFO] Command: {' '.join(cmd)}")
    print("-" * 60)

    try:
        subprocess.run(cmd, check=True)
        print(f"[SUCCESS] Sketch created: {sketch_file}.msh")
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Mash sketch failed with exit code {e.returncode}")
        sys.exit(1)
    except Exception as e:
        print(f"[ERROR] An unexpected error occurred: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()