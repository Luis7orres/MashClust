#!/usr/bin/env python3
"""
Script 3: Cluster Genomes and Select Representatives (Multi-Reference Support)
Low-RAM: Uses disk-based sparse pair storage to handle large distance matrices.
"""
import sys
import argparse
import json
import numpy as np
import re
import os
import tempfile
from pathlib import Path
from collections import defaultdict

def extract_genome_id(path_str):
    """Extract clean ID (GCF/GCA) or fallback to folder name"""
    path = Path(path_str)
    match = re.search(r'(GC[FA]_\d{9}[._]\d+)', path.name)
    if match: return match.group(1)
    match = re.search(r'(GC[FA]_\d{9}[._]\d+)', path.parent.name)
    if match: return match.group(1)
    return path.parent.name

def normalize_id(acc_id):
    """Normalize Accession ID by replacing dots with underscores."""
    return acc_id.replace('.', '_')

def parse_matrix_to_pairs(dist_file, distance_threshold, pairs_file):
    """
    PASS 1: Stream the distance matrix line by line.
    Write only pairs BELOW the threshold to a temp file on disk.
    Also collect genome_paths mapping. Never loads the full matrix into RAM.
    Returns: genome_paths dict, total genomes count.
    """
    genome_paths = {}
    genome_ids_ordered = []
    pairs_written = 0

    print("[INFO] Pass 1/2 — Streaming distance matrix, writing sparse pairs to disk...")

    with open(dist_file, 'r') as f_in, open(pairs_file, 'w') as f_out:
        # Read header to build genome list
        header_line = f_in.readline().strip()
        header = header_line.split('\t')[1:]
        all_genomes_list = [p.strip() for p in header]

        for p in all_genomes_list:
            gid = extract_genome_id(p)
            genome_paths[gid] = p
        genome_ids_ordered = list(genome_paths.keys())

        # Stream rows
        for i, line in enumerate(f_in):
            if i % 5000 == 0 and i > 0:
                print(f"    [INFO] Processed {i:,} / {len(genome_ids_ordered):,} rows, "
                      f"{pairs_written:,} pairs written...", flush=True)

            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                continue

            query_id = extract_genome_id(parts[0])
            dists = parts[1:]

            for j, d_str in enumerate(dists):
                try:
                    d = float(d_str)
                    target_id = genome_ids_ordered[j]
                    # Only write pairs that are similar enough AND avoid self-pairs
                    # Write only j > i to avoid duplicates (upper triangle)
                    if query_id != target_id and d <= distance_threshold and j > i:
                        f_out.write(f"{query_id}\t{target_id}\n")
                        pairs_written += 1
                except (ValueError, IndexError):
                    continue

    print(f"[INFO] Pass 1 complete: {len(genome_ids_ordered):,} genomes, "
          f"{pairs_written:,} pairs below threshold.")
    return genome_paths, genome_ids_ordered

def build_neighbors_from_pairs(pairs_file):
    """
    PASS 2: Read sparse pairs from disk and build the neighbors dict.
    """
    print("[INFO] Pass 2/2 — Building neighbor graph from sparse pairs...")
    neighbors = defaultdict(set)
    count = 0

    with open(pairs_file, 'r') as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) == 2:
                a, b = parts[0], parts[1]
                neighbors[a].add(b)
                neighbors[b].add(a)
                count += 1

    print(f"[INFO] Pass 2 complete: {count:,} pairs loaded, "
          f"{len(neighbors):,} genomes with at least one neighbor.")
    return neighbors

def main():
    parser = argparse.ArgumentParser(description="Step 3: Cluster Genomes (Low-RAM Version)")
    parser.add_argument("input_dir", help="Directory with distances.txt or path to distances.txt file")
    parser.add_argument("-o", "--output", required=True, help="Output file (representatives.txt)")
    parser.add_argument("-t", "--threshold", type=float, default=0.9997)
    parser.add_argument("-n", "--num-representatives", type=int, default=5)
    parser.add_argument("--reference", nargs='+', help="Reference Accession(s) to protect")
    parser.add_argument("--no-reference-protection", action="store_true")
    args = parser.parse_args()

    # Output setup — identical to original
    out_file = Path(args.output)
    out_dir = out_file.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    # Locate distances.txt — identical to original
    input_path = Path(args.input_dir)
    dist_file = input_path if input_path.is_file() else input_path / "distances.txt"

    if not dist_file.exists():
        print(f"[ERROR] Distance file not found at: {dist_file}")
        sys.exit(1)

    print(f"[INFO] Threshold: {args.threshold} | Max Reps: {args.num_representatives}")
    distance_threshold = 1.0 - args.threshold

    # ----------------------------------------------------------------
    # Two-pass streaming
    # ----------------------------------------------------------------
    # Temp file for sparse pairs — written in the output dir to use Store space
    pairs_file = out_dir / ".sparse_pairs_tmp.tsv"

    try:
        # PASS 1: Stream matrix → write sparse pairs to disk
        genome_paths, genome_ids_ordered = parse_matrix_to_pairs(
            dist_file, distance_threshold, pairs_file
        )

        # PASS 2: Load only the sparse pairs into memory
        neighbors = build_neighbors_from_pairs(pairs_file)

    finally:
        # Always clean up temp file
        if pairs_file.exists():
            pairs_file.unlink()
            print("[INFO] Temporary pairs file cleaned up.")

    # 2. Identify References
    ref_ids = set()
    if not args.no_reference_protection and args.reference:
        for ref_input in args.reference:
            ref_input_clean = normalize_id(ref_input)
            found = False
            for gid in genome_paths.keys():
                if ref_input_clean in normalize_id(gid):
                    ref_ids.add(gid)
                    print(f"[INFO] Reference genome identified: {gid} (matched {ref_input})")
                    found = True
                    break
            if not found:
                print(f"[WARNING] Reference {ref_input} not found in dataset.")

    # 3. Greedy Clustering
    print("[INFO] Running Greedy Clustering...")
    clusters = []
    assigned = set()

    genome_connectivity = [(g, len(neighbors[g])) for g in genome_paths.keys()]

    if ref_ids:
        refs = [x for x in genome_connectivity if x[0] in ref_ids]
        others = [x for x in genome_connectivity if x[0] not in ref_ids]
        refs.sort(key=lambda x: x[1], reverse=True)
        others.sort(key=lambda x: x[1], reverse=True)
        genome_connectivity = refs + others
        print(f"[INFO] {len(refs)} references prioritized in clustering queue.")
    else:
        genome_connectivity.sort(key=lambda x: x[1], reverse=True)

    for g, _ in genome_connectivity:
        if g in assigned:
            continue

        cluster = {g}
        assigned.add(g)

        for neighbor in neighbors[g]:
            if neighbor not in assigned:
                cluster.add(neighbor)
                assigned.add(neighbor)

        clusters.append(cluster)

    print(f"[INFO] Created {len(clusters)} clusters.")

    # 4. Select Representatives
    representatives = []
    np.random.seed(42)

    for cluster in clusters:
        cluster_list = list(cluster)

        if len(cluster_list) <= args.num_representatives:
            representatives.extend(cluster_list)
        else:
            selected = []
            pool = list(cluster_list)

            if ref_ids:
                refs_in_cluster = [m for m in pool if m in ref_ids]
                for r in refs_in_cluster:
                    selected.append(r)
                    pool.remove(r)
                    if len(selected) >= args.num_representatives:
                        break

            slots_left = args.num_representatives - len(selected)
            if slots_left > 0:
                picked = np.random.choice(pool, size=slots_left, replace=False)
                selected.extend(picked)

            representatives.extend(selected)

    # 5. Output
    # A) Representatives list
    with open(out_file, "w") as f:
        for r in representatives:
            f.write(f"{genome_paths[r]}\n")

    # B) JSON data — same structure as original, compatible with scripts 4 and 5
    json_path = out_dir / "clustering_data.json"
    json_data = {
        'clusters': [list(c) for c in clusters],
        'neighbors': {k: list(v) for k, v in neighbors.items()},
        'genome_names': genome_paths,
        'identity_threshold': args.threshold
    }

    with open(json_path, "w") as f:
        json.dump(json_data, f)

    print(f"[SUCCESS] Selected {len(representatives)} representatives.")
    print(f"[SUCCESS] Representatives written to: {out_file}")
    print(f"[SUCCESS] Clustering data written to: {json_path}")

if __name__ == "__main__":
    main()