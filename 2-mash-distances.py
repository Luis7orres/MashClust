#!/usr/bin/env python3
"""
Script 2: Calculate distances, perform clustering, and select representatives
Reads sketch file, calculates distances, clusters, and selects representatives
"""

import subprocess
import sys
from pathlib import Path
import numpy as np
from collections import defaultdict
import argparse
import re
from datetime import datetime
import json


class Logger:
    """Logger that writes to both stdout and file"""
    def __init__(self, log_file):
        self.terminal = sys.stdout
        self.log = open(log_file, 'w', buffering=1)
        self.log.write(f"=== Distance log started at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} ===\n\n")
    
    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
    
    def flush(self):
        self.terminal.flush()
        self.log.flush()
    
    def close(self):
        self.log.close()


def extract_genome_id(path_str):
    """Extract genome ID from path"""
    path = Path(path_str)
    filename = path.name
    
    match = re.search(r'(GC[FA]_\d{9}_\d+)', filename)
    if match:
        return match.group(1)
    
    match = re.search(r'(GC[FA]_\d{9}_\d+)', path.parent.name)
    if match:
        return match.group(1)
    
    return str(path)


def calculate_distances(sketch_file, output_file, threads=1):
    """Calculate all-vs-all distance matrix"""
    
    print(f"\n{'='*70}")
    print(f"CALCULATING MASH DISTANCES")
    print(f"{'='*70}")
    print(f"[INFO] Sketch file: {sketch_file}.msh")
    print(f"[INFO] Output file: {output_file}")
    print(f"[INFO] Threads: {threads}")
    
    if not Path(f"{sketch_file}.msh").exists():
        print(f"[ERROR] Sketch file not found: {sketch_file}.msh")
        print(f"[ERROR] Run mash_sketch.py first!")
        sys.exit(1)
    
    cmd = [
        "mash", "dist",
        "-t",  # Tabular output (matrix format)
        "-p", str(threads),
        f"{sketch_file}.msh",
        f"{sketch_file}.msh"
    ]
    
    print(f"\n[INFO] Running Mash dist...")
    print(f"[INFO] This may take 10-20 minutes...")
    print(f"[CMD] {' '.join(cmd)}\n")
    
    try:
        with open(output_file, 'w') as f:
            result = subprocess.run(cmd, stdout=f, check=True, stderr=subprocess.PIPE, text=True)
        
        # Show stderr for progress
        if result.stderr:
            print(result.stderr)
        
        # Check output file
        output_path = Path(output_file)
        if output_path.exists():
            size_mb = output_path.stat().st_size / (1024 * 1024)
            print(f"\n[SUCCESS] Distances calculated: {output_file}")
            print(f"[INFO] Distance file size: {size_mb:.2f} MB")
            
            # Count lines
            with open(output_file, 'r') as f:
                line_count = sum(1 for _ in f)
            print(f"[INFO] Matrix lines: {line_count}")
        
    except subprocess.CalledProcessError as e:
        print(f"\n[ERROR] Mash dist failed!")
        print(f"[ERROR] {e.stderr}")
        sys.exit(1)


def parse_mash_distances(dist_file, identity_threshold=0.998):
    """Parse Mash distance matrix and build neighbor graph"""
    
    print(f"\n{'='*70}")
    print(f"PARSING DISTANCE MATRIX")
    print(f"{'='*70}")
    print(f"[INFO] Identity threshold: {identity_threshold*100:.2f}%")
    print(f"[INFO] Distance threshold: {1-identity_threshold:.6f}")
    
    distance_threshold = 1 - identity_threshold
    neighbors = defaultdict(set)
    genome_names = {}
    
    line_count = 0
    similar_count = 0
    min_distance = float('inf')
    max_distance = 0
    
    with open(dist_file, 'r') as f:
        # Read header
        header_line = f.readline().strip()
        
        if not header_line.startswith('#query'):
            print(f"[ERROR] Invalid distance file format!")
            print(f"[ERROR] Expected header starting with '#query'")
            sys.exit(1)
        
        header_parts = header_line.split('\t')
        genome_order = header_parts[1:]
        
        print(f"[INFO] Matrix dimensions: {len(genome_order)} x {len(genome_order)}")
        print(f"[INFO] First 3 genomes:")
        for i, g in enumerate(genome_order[:3], 1):
            genome_id = extract_genome_id(g)
            print(f"  {i}. {genome_id}")
            genome_names[genome_id] = Path(g)
        
        # Store all genome names
        for genome_path in genome_order:
            genome_id = extract_genome_id(genome_path)
            genome_names[genome_id] = Path(genome_path)
        
        print(f"[INFO] Total genomes: {len(genome_names)}")
        print(f"\n[INFO] Parsing distance matrix...")
        
        # Parse distance rows
        for row_idx, line in enumerate(f):
            if not line.strip():
                continue
            
            parts = line.strip().split('\t')
            
            if len(parts) < 2:
                continue
            
            query_path = parts[0]
            query_id = extract_genome_id(query_path)
            
            if query_id not in genome_names:
                genome_names[query_id] = Path(query_path)
            
            distances = parts[1:]
            
            # Process distances
            for col_idx, dist_str in enumerate(distances[:len(genome_order)]):
                try:
                    distance = float(dist_str)
                except ValueError:
                    continue
                
                target_path = genome_order[col_idx]
                target_id = extract_genome_id(target_path)
                
                # Track distance statistics
                if query_id != target_id:
                    min_distance = min(min_distance, distance)
                    max_distance = max(max_distance, distance)
                
                # Build neighbor graph
                if distance <= distance_threshold and query_id != target_id:
                    neighbors[query_id].add(target_id)
                    neighbors[target_id].add(query_id)
                    similar_count += 1
            
            line_count += 1
            
            if line_count % 500 == 0:
                print(f"[PROGRESS] Processed {line_count} rows...")
    
    similar_count = similar_count // 2  # Each pair counted twice
    
    print(f"\n[INFO] Parsing complete!")
    print(f"[INFO] Rows processed: {line_count}")
    
    if min_distance != float('inf'):
        print(f"[INFO] Distance range: {min_distance:.6f} to {max_distance:.6f}")
        print(f"[INFO] Identity range: {(1-max_distance)*100:.3f}% to {(1-min_distance)*100:.3f}%")
    
    print(f"[INFO] Similar pairs (≤{distance_threshold:.6f}): {similar_count}")
    print(f"[INFO] Genomes with neighbors: {len([g for g in genome_names if g in neighbors])}")
    print(f"[INFO] Isolated genomes: {len([g for g in genome_names if g not in neighbors])}")
    
    if similar_count == 0 and len(genome_names) > 1:
        print(f"\n[WARNING] No similar genome pairs found!")
        print(f"[WARNING] Current threshold may be too strict: {identity_threshold*100:.2f}%")
        if min_distance != float('inf'):
            suggested_threshold = 1 - (min_distance * 1.1)  # 10% margin
            print(f"[SUGGESTION] Try a lower threshold like: {suggested_threshold:.4f}")
    
    return neighbors, genome_names


def greedy_clustering(neighbors, all_genomes):
    """Perform greedy clustering"""
    
    print(f"\n{'='*70}")
    print(f"CLUSTERING GENOMES")
    print(f"{'='*70}")
    
    clusters = []
    assigned = set()
    
    all_genomes = list(all_genomes)
    
    # Sort by connectivity (most connected first)
    genome_connectivity = [(g, len(neighbors.get(g, set()))) for g in all_genomes]
    genome_connectivity.sort(key=lambda x: x[1], reverse=True)
    
    print(f"[INFO] Starting greedy clustering...")
    
    for genome, conn in genome_connectivity:
        if genome in assigned:
            continue
        
        # Start new cluster
        cluster = {genome}
        assigned.add(genome)
        
        # Add all unassigned neighbors
        for neighbor in neighbors.get(genome, set()):
            if neighbor not in assigned:
                cluster.add(neighbor)
                assigned.add(neighbor)
        
        clusters.append(cluster)
    
    print(f"[INFO] Clustering complete!")
    print(f"[INFO] Total clusters: {len(clusters)}")
    
    # Statistics
    cluster_sizes = [len(c) for c in clusters]
    print(f"\nCluster statistics:")
    print(f"  Mean size:     {np.mean(cluster_sizes):.2f}")
    print(f"  Median size:   {np.median(cluster_sizes):.0f}")
    print(f"  Largest:       {max(cluster_sizes)}")
    print(f"  Smallest:      {min(cluster_sizes)}")
    print(f"  Singletons:    {sum(1 for s in cluster_sizes if s == 1)}")
    
    # Show top clusters
    clusters_sorted = sorted(clusters, key=len, reverse=True)
    print(f"\nTop 10 largest clusters:")
    for i, cluster in enumerate(clusters_sorted[:10], 1):
        sample_genomes = list(cluster)[:3]
        sample_str = ', '.join(str(g)[:40] for g in sample_genomes)
        if len(cluster) > 3:
            sample_str += ', ...'
        print(f"  {i:2d}. Size {len(cluster):5d}: {sample_str}")
    
    return clusters


def select_representatives(clusters, max_per_cluster=5, seed=42):
    """Select representative genomes from each cluster"""
    
    print(f"\n{'='*70}")
    print(f"SELECTING REPRESENTATIVES")
    print(f"{'='*70}")
    print(f"[INFO] Max representatives per cluster: {max_per_cluster}")
    
    np.random.seed(seed)
    representatives = []
    
    large_clusters = 0
    
    for i, cluster in enumerate(clusters):
        cluster_list = list(cluster)
        
        if len(cluster_list) <= max_per_cluster:
            selected = cluster_list
        else:
            selected = np.random.choice(
                cluster_list, 
                size=max_per_cluster, 
                replace=False
            ).tolist()
            large_clusters += 1
        
        representatives.extend(selected)
    
    print(f"[INFO] Representatives selected: {len(representatives)}")
    print(f"[INFO] Large clusters (>{max_per_cluster}): {large_clusters}")
    print(f"[INFO] Original genomes: {sum(len(c) for c in clusters)}")
    print(f"[INFO] Reduction: {(1 - len(representatives)/sum(len(c) for c in clusters))*100:.1f}%")
    
    return representatives


def write_outputs(clusters, representatives, genome_names, output_dir):
    """Write all output files"""
    
    print(f"\n{'='*70}")
    print(f"WRITING OUTPUT FILES")
    print(f"{'='*70}")
    
    output_dir = Path(output_dir)
    
    # 1. Representatives list
    rep_file = output_dir / "representatives.txt"
    missing = []
    
    with open(rep_file, 'w') as f:
        for rep in representatives:
            if rep in genome_names:
                f.write(f"{genome_names[rep]}\n")
            else:
                missing.append(rep)
    
    if missing:
        print(f"[WARNING] {len(missing)} representatives not found in genome paths")
    
    print(f"[INFO] Representatives: {rep_file}")
    print(f"        ({len(representatives)} genomes)")
    
    # 2. Clusters file
    clusters_file = output_dir / "clusters.txt"
    
    with open(clusters_file, 'w') as f:
        f.write("Cluster_ID\tSize\tMembers\n")
        for i, cluster in enumerate(sorted(clusters, key=len, reverse=True), 1):
            members = ','.join(str(g) for g in sorted(cluster))
            f.write(f"Cluster_{i}\t{len(cluster)}\t{members}\n")
    
    print(f"[INFO] Clusters:        {clusters_file}")
    print(f"        ({len(clusters)} clusters)")
    
    # 3. Cluster summary
    summary_file = output_dir / "cluster_summary.txt"
    
    cluster_sizes = [len(c) for c in clusters]
    
    with open(summary_file, 'w') as f:
        f.write("CLUSTERING SUMMARY\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Total genomes:        {sum(cluster_sizes)}\n")
        f.write(f"Total clusters:       {len(clusters)}\n")
        f.write(f"Representatives:      {len(representatives)}\n")
        f.write(f"Reduction:            {(1-len(representatives)/sum(cluster_sizes))*100:.1f}%\n\n")
        f.write(f"Cluster size statistics:\n")
        f.write(f"  Mean:               {np.mean(cluster_sizes):.2f}\n")
        f.write(f"  Median:             {np.median(cluster_sizes):.0f}\n")
        f.write(f"  Min:                {min(cluster_sizes)}\n")
        f.write(f"  Max:                {max(cluster_sizes)}\n")
        f.write(f"  Singletons:         {sum(1 for s in cluster_sizes if s == 1)}\n")
    
    print(f"[INFO] Summary:         {summary_file}")
    
    return rep_file, clusters_file, summary_file


def save_clustering_data(clusters, neighbors, genome_names, output_dir, identity_threshold):
    """Save clustering data in JSON for visualization"""
    
    print(f"[INFO] Saving data for visualization...")
    
    # Convert to serializable format
    clusters_list = [list(c) for c in clusters]
    neighbors_dict = {k: list(v) for k, v in neighbors.items()}
    genome_names_dict = {k: str(v) for k, v in genome_names.items()}
    
    data = {
        'clusters': clusters_list,
        'neighbors': neighbors_dict,
        'genome_names': genome_names_dict,
        'identity_threshold': identity_threshold,
        'n_clusters': len(clusters),
        'n_genomes': len(genome_names),
        'timestamp': datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    }
    
    output_file = output_dir / 'clustering_data.json'
    with open(output_file, 'w') as f:
        json.dump(data, f, indent=2)
    
    print(f"[INFO] Visualization data: {output_file}")
    
    return output_file


def main():
    parser = argparse.ArgumentParser(
        description="Step 2: Calculate distances and perform clustering",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage (uses sketch from mash_output/)
  python mash_distances.py mash_output
  
  # Custom threshold
  python mash_distances.py mash_output -t 0.999 -n 10
  
  # With more threads
  python mash_distances.py mash_output --threads 8

Next step:
  python mash_visualize.py mash_output
        """
    )
    
    parser.add_argument("output_dir", help="Directory containing sketch file")
    parser.add_argument("-t", "--threshold", type=float, default=0.9995,
                       help="Identity threshold (default: 0.9995 = 99.95%%)")
    parser.add_argument("-n", "--num-representatives", type=int, default=5,
                       help="Representatives per cluster (default: 5)")
    parser.add_argument("--threads", type=int, default=1,
                       help="CPU threads (default: 1, use -1 for all)")
    parser.add_argument("--seed", type=int, default=42,
                       help="Random seed for reproducibility")
    
    args = parser.parse_args()
    
    # Handle threads
    threads = args.threads
    if threads == -1:
        import multiprocessing
        threads = multiprocessing.cpu_count()
        print(f"[INFO] Using all {threads} CPU threads")
    
    output_dir = Path(args.output_dir)
    
    if not output_dir.exists():
        print(f"[ERROR] Directory not found: {output_dir}")
        print(f"[ERROR] Run mash_sketch.py first!")
        sys.exit(1)
    
    # Setup logging
    log_file = output_dir / f"distances_log_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
    logger = Logger(log_file)
    sys.stdout = logger
    sys.stderr = logger
    
    print(f"[INFO] Log file: {log_file}\n")
    
    sketch_file = output_dir / "genomes_sketch"
    dist_file = output_dir / "distances.txt"
    
    try:
        print("="*70)
        print("MASH DISTANCES AND CLUSTERING - STEP 2/3")
        print("="*70)
        
        # Calculate distances
        calculate_distances(sketch_file, dist_file, threads)
        
        # Parse distances
        neighbors, genome_names = parse_mash_distances(dist_file, args.threshold)
        
        if not genome_names:
            print("[ERROR] No genomes found!")
            sys.exit(1)
        
        # Cluster
        clusters = greedy_clustering(neighbors, genome_names.keys())
        
        if not clusters:
            print("[ERROR] No clusters created!")
            sys.exit(1)
        
        # Select representatives
        representatives = select_representatives(clusters, args.num_representatives, args.seed)
        
        # Write outputs
        rep_file, clusters_file, summary_file = write_outputs(
            clusters, representatives, genome_names, output_dir
        )
        
        # Save for visualization
        data_file = save_clustering_data(
            clusters, neighbors, genome_names, output_dir, args.threshold
        )
        
        # Final summary
        print(f"\n{'='*70}")
        print("SUMMARY")
        print(f"{'='*70}")
        print(f"Original genomes:       {len(genome_names)}")
        print(f"Clusters:               {len(clusters)}")
        print(f"Representatives:        {len(representatives)}")
        print(f"Reduction:              {(1-len(representatives)/len(genome_names))*100:.1f}%")
        print(f"Threshold:              {args.threshold*100:.2f}% identity")
        print(f"\nOutput files:")
        print(f"  - {rep_file.name}")
        print(f"  - {clusters_file.name}")
        print(f"  - {summary_file.name}")
        print(f"  - {dist_file.name}")
        print(f"  - {data_file.name}")
        print(f"  - {log_file.name}")
        print(f"\nNext step:")
        print(f"  python mash_visualize.py {args.output_dir}")
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