#!/usr/bin/env python3
"""
Script 1: Create Mash sketch from genome files
Finds genome files and creates a Mash sketch
"""

import subprocess
import os
import sys
from pathlib import Path
import argparse
import tempfile
import re
from datetime import datetime


class Logger:
    """Logger that writes to both stdout and file"""
    def __init__(self, log_file):
        self.terminal = sys.stdout
        self.log = open(log_file, 'w', buffering=1)
        self.log.write(f"=== Sketch log started at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} ===\n\n")
    
    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
    
    def flush(self):
        self.terminal.flush()
        self.log.flush()
    
    def close(self):
        self.log.close()


def find_genome_files(base_dir, file_pattern="GENOME_*.fna", parent_filter="mycobacterium_tuberculosis"):
    """Find genome files in NCBI directory structure"""
    base_path = Path(base_dir)
    genome_files = []
    filtered_out = []
    
    print(f"[INFO] Searching for genome files matching '{file_pattern}'")
    print(f"[INFO] FILTER: Only in directories starting with '{parent_filter}*'")
    
    # Search recursively
    all_matches = []
    for pattern in [file_pattern, file_pattern + ".gz"]:
        all_matches.extend(base_path.rglob(pattern))
    
    # Filter by parent directory name
    for genome_file in all_matches:
        parent_dir = genome_file.parent.name.lower()
        
        if parent_dir.startswith(parent_filter.lower()):
            genome_files.append(genome_file)
        else:
            filtered_out.append(genome_file)
    
    genome_files = sorted(set(genome_files))
    
    print(f"[INFO] Found {len(genome_files)} M. tuberculosis genome files")
    print(f"[INFO] Filtered out {len(filtered_out)} genomes from other species")
    
    if not genome_files:
        print(f"[ERROR] No genome files found matching filter '{parent_filter}*'")
        sys.exit(1)
    
    # Show samples
    if len(genome_files) > 0:
        print(f"[INFO] Sample of M. tuberculosis genomes:")
        for f in genome_files[:5]:
            print(f"       {f.parent.name}/{f.name}")
        if len(genome_files) > 5:
            print(f"       ... and {len(genome_files) - 5} more")
    
    return genome_files


def save_genome_list(genome_files, output_file):
    """Save list of genome files for later use"""
    print(f"[INFO] Saving genome file list...")
    
    with open(output_file, 'w') as f:
        for genome_file in genome_files:
            f.write(f"{genome_file}\n")
    
    print(f"[INFO] Saved {len(genome_files)} genome paths to: {output_file}")


def create_mash_sketch(genome_files, output_sketch, sketch_size=100000, kmer_size=31, threads=1):
    """Create Mash sketch for all genomes"""
    
    print(f"\n{'='*70}")
    print(f"CREATING MASH SKETCH")
    print(f"{'='*70}")
    print(f"[INFO] Number of genomes: {len(genome_files)}")
    print(f"[INFO] Parameters:")
    print(f"       k-mer size:    {kmer_size}")
    print(f"       sketch size:   {sketch_size}")
    print(f"       threads:       {threads}")
    print(f"[INFO] Output: {output_sketch}.msh")
    
    # Create temporary file list
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as tmp:
        tmp_file = tmp.name
        for f in genome_files:
            tmp.write(f"{f}\n")
    
    try:
        cmd = [
            "mash", "sketch",
            "-s", str(sketch_size),
            "-k", str(kmer_size),
            "-p", str(threads),
            "-o", str(output_sketch),
            "-l", tmp_file
        ]
        
        print(f"\n[INFO] Running Mash sketch...")
        print(f"[INFO] This may take 15-30 minutes depending on dataset size...")
        print(f"[CMD] {' '.join(cmd)}\n")
        
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        # Show Mash output
        if result.stderr:
            print(result.stderr)
        
        print(f"\n[SUCCESS] Sketch created: {output_sketch}.msh")
        
        # Check file size
        sketch_file = Path(f"{output_sketch}.msh")
        if sketch_file.exists():
            size_mb = sketch_file.stat().st_size / (1024 * 1024)
            print(f"[INFO] Sketch file size: {size_mb:.2f} MB")
        
        return len(genome_files)
        
    except subprocess.CalledProcessError as e:
        print(f"\n[ERROR] Mash sketch failed!")
        print(f"[ERROR] {e.stderr}")
        sys.exit(1)
    finally:
        if os.path.exists(tmp_file):
            os.remove(tmp_file)


def main():
    parser = argparse.ArgumentParser(
        description="Step 1: Create Mash sketch from genome files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  python mash_sketch.py /data/genomes -o mash_output
  
  # Custom parameters
  python mash_sketch.py /data/genomes -o mash_output -s 50000 -k 21 --threads 8
  
  # Custom file pattern
  python mash_sketch.py /data/genomes -o output -p "*.fasta" -f "tuberculosis"

Next step:
  python mash_distances.py mash_output
        """
    )
    
    parser.add_argument("genome_dir", help="Base directory containing genome subdirectories")
    parser.add_argument("-o", "--output-dir", default="mash_output", 
                       help="Output directory (default: mash_output)")
    parser.add_argument("-s", "--sketch-size", type=int, default=100000,
                       help="Mash sketch size (default: 100000)")
    parser.add_argument("-k", "--kmer-size", type=int, default=31,
                       help="K-mer size (default: 31)")
    parser.add_argument("-p", "--pattern", default="GENOME_*.fna",
                       help="File pattern to search (default: GENOME_*.fna)")
    parser.add_argument("-f", "--filter", default="mycobacterium_tuberculosis",
                       help="Parent directory filter (default: mycobacterium_tuberculosis)")
    parser.add_argument("--threads", type=int, default=1, 
                       help="CPU threads to use (default: 1, use -1 for all)")
    
    args = parser.parse_args()
    
    # Handle threads
    threads = args.threads
    if threads == -1:
        import multiprocessing
        threads = multiprocessing.cpu_count()
        print(f"[INFO] Using all {threads} CPU threads")
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Setup logging
    log_file = output_dir / f"sketch_log_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
    logger = Logger(log_file)
    sys.stdout = logger
    sys.stderr = logger
    
    print(f"[INFO] Log file: {log_file}\n")
    
    sketch_file = output_dir / "genomes_sketch"
    genome_list_file = output_dir / "genome_list.txt"
    
    try:
        print("="*70)
        print("MASH SKETCH CREATION - STEP 1/3")
        print("="*70)
        
        # Find genome files
        genome_files = find_genome_files(args.genome_dir, args.pattern, args.filter)
        
        # Save genome list for reference
        save_genome_list(genome_files, genome_list_file)
        
        # Create sketch
        num_genomes = create_mash_sketch(
            genome_files, 
            sketch_file, 
            args.sketch_size, 
            args.kmer_size, 
            threads
        )
        
        # Summary
        print(f"\n{'='*70}")
        print("SUMMARY")
        print(f"{'='*70}")
        print(f"Genomes processed:  {num_genomes}")
        print(f"Sketch size:        {args.sketch_size}")
        print(f"K-mer size:         {args.kmer_size}")
        print(f"Threads used:       {threads}")
        print(f"\nOutput files:")
        print(f"  - Sketch:       {sketch_file}.msh")
        print(f"  - Genome list:  {genome_list_file}")
        print(f"  - Log file:     {log_file}")
        print(f"\nNext step:")
        print(f"  python 2-mash-distances.py {args.output_dir}")
        print("="*70)
        
    except Exception as e:
        print(f"\n[ERROR] Failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    finally:
        logger.close()
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__


if __name__ == "__main__":
    main()