#!/usr/bin/env python3
"""
Script 3: Cluster Genomes and Select Representatives (GCF Priority Version)
"""
import sys
import argparse
import json
import numpy as np
import re
from pathlib import Path
from collections import defaultdict

def extract_genome_id(path_str):
    """Extract clean ID (GCF/GCA) or fallback to folder name"""
    path = Path(path_str)
    # Attempt 1: File name
    match = re.search(r'(GC[FA]_\d{9}_\d+)', path.name)
    if match: return match.group(1)
    # Attempt 2: Parent folder name
    match = re.search(r'(GC[FA]_\d{9}_\d+)', path.parent.name)
    if match: return match.group(1)
    # Fallback
    return path.parent.name

def main():
    parser = argparse.ArgumentParser(description="Step 3: Cluster Genomes")
    parser.add_argument("input_dir", help="Directory with distances.txt")
    parser.add_argument("-o", "--output-dir", required=True)
    parser.add_argument("-t", "--threshold", type=float, default=0.9997)
    parser.add_argument("-n", "--num-representatives", type=int, default=5)
    parser.add_argument("--reference", type=str, help="Reference Accession to protect")
    parser.add_argument("--no-reference-protection", action="store_true")
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    dist_file = Path(args.input_dir) / "distances.txt"

    print(f"[INFO] Threshold: {args.threshold} | Max Reps: {args.num_representatives}")
    
    # 1. Parse Distance Matrix
    neighbors = defaultdict(set)
    genome_paths = {} # Map ID -> Full Path
    
    distance_threshold = 1.0 - args.threshold
    
    print("[INFO] Parsing distance matrix...")
    with open(dist_file, 'r') as f:
        # Read header
        header = f.readline().strip().split('\t')[1:]
        all_genomes_list = [p.strip() for p in header]
        
        # Map IDs to paths
        for p in all_genomes_list:
            gid = extract_genome_id(p)
            genome_paths[gid] = p
            
        genome_ids_ordered = list(genome_paths.keys())

        # Read rows
        for i, line in enumerate(f):
            parts = line.strip().split('\t')
            if len(parts) < 2: continue
            
            query_id = extract_genome_id(parts[0])
            dists = parts[1:]
            
            for j, d_str in enumerate(dists):
                try:
                    d = float(d_str)
                    target_id = genome_ids_ordered[j]
                    
                    if query_id != target_id and d <= distance_threshold:
                        neighbors[query_id].add(target_id)
                        neighbors[target_id].add(query_id)
                except ValueError:
                    continue

    # 2. Identify Reference
    ref_id = None
    if not args.no_reference_protection and args.reference:
        for gid in genome_paths.keys():
            if args.reference in gid:
                ref_id = gid
                print(f"[INFO] Reference genome identified: {ref_id}")
                break
        if not ref_id:
            print(f"[WARNING] Reference {args.reference} not found in dataset.")

    # 3. Greedy Clustering
    print("[INFO] Running Greedy Clustering...")
    clusters = []
    assigned = set()
    
    # Sort genomes: Reference first (forced), then by connectivity (most neighbors first)
    genome_connectivity = [(g, len(neighbors[g])) for g in genome_paths.keys()]
    
    if ref_id:
        genome_connectivity = [x for x in genome_connectivity if x[0] != ref_id]
        genome_connectivity.insert(0, (ref_id, len(neighbors[ref_id])))
        print("[INFO] Reference prioritized in clustering queue.")
    else:
        genome_connectivity.sort(key=lambda x: x[1], reverse=True)

    for g, _ in genome_connectivity:
        if g in assigned:
            continue
            
        # New cluster
        cluster = {g}
        assigned.add(g)
        
        # Add neighbors
        for neighbor in neighbors[g]:
            if neighbor not in assigned:
                cluster.add(neighbor)
                assigned.add(neighbor)
        
        clusters.append(cluster)

    print(f"[INFO] Created {len(clusters)} clusters.")

    # 4. Select Representatives (Prioritizing GCF)
    representatives = []
    np.random.seed(42) # Reproducibility
    
    for cluster in clusters:
        cluster_list = list(cluster)
        
        # If cluster size is smaller than required reps, take all
        if len(cluster_list) <= args.num_representatives:
            representatives.extend(cluster_list)
        else:
            selected = []
            pool = list(cluster_list)
            
            # A) Pick reference first if present (Reference protection)
            if ref_id and ref_id in pool:
                selected.append(ref_id)
                pool.remove(ref_id)
            
            # B) Prioritize GCFs in the remaining pool
            # Separate pool into GCF and GCA groups
            gcfs = [g for g in pool if g.startswith('GCF')]
            gcas = [g for g in pool if not g.startswith('GCF')]
            
            # Shuffle each group internally to maintain randomness among equals
            np.random.shuffle(gcfs)
            np.random.shuffle(gcas)
            
            # Reconstruct ordered pool: GCFs always come before GCAs
            ordered_pool = gcfs + gcas
            
            # C) Fill remaining slots following priority order (GCF > GCA)
            slots_left = args.num_representatives - len(selected)
            if slots_left > 0:
                # Pick the first N from the prioritized pool
                picked = ordered_pool[:slots_left]
                selected.extend(picked)
            
            representatives.extend(selected)

    # 5. Output
    # A) Representatives list (full paths)
    rep_file = out_dir / "representatives.txt"
    with open(rep_file, "w") as f:
        for r in representatives:
            f.write(f"{genome_paths[r]}\n")
            
    # B) JSON data for visualization
    json_data = {
        'clusters': [list(c) for c in clusters],
        'neighbors': {k: list(v) for k, v in neighbors.items()},
        'genome_names': genome_paths,
        'identity_threshold': args.threshold
    }
    
    with open(out_dir / "clustering_data.json", "w") as f:
        json.dump(json_data, f)
        
    print(f"[SUCCESS] Selected {len(representatives)} representatives (Prioritizing GCF).")

if __name__ == "__main__":
    main()