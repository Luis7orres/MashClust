#!/usr/bin/env python
"""
Script 5: Finalize Results (The "Perfect Integration" approach)
Uses dataset-manager.py directly to guarantee 100% compatibility with Pipeline 1.
"""
import argparse
import shutil
import re
import os
import sys
import subprocess
import time
from pathlib import Path

def normalize_accession(path_str):
    """Extracts GCA_000000000.1 from any string."""
    pattern = r'(GC[FA])_(\d{9})[_.](\d+)'
    match = re.search(pattern, path_str)
    if match:
        return f"{match.group(1)}_{match.group(2)}.{match.group(3)}"
    return None

def main():
    parser = argparse.ArgumentParser(description="Step 5: Finalize Results")
    parser.add_argument("--non-targets", required=True)
    parser.add_argument("--representatives", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--acc-file", required=True, help="Path for selected_accessions.txt")
    parser.add_argument("--genomes-subdir", default="uncompressed")
    parser.add_argument("--dataset-manager", required=True, help="Absolute path to dataset-manager.py")
    args = parser.parse_args()

    print("\n" + "="*60)
    print(" STARTING INTEGRATION WITH DATASET-MANAGER")
    print("="*60)

    # 1. Verify tools
    datasets_exe = shutil.which("datasets")
    if not datasets_exe:
        print("[CRITICAL ERROR] 'datasets' binary not found.")
        return

    dm_script = Path(args.dataset_manager)
    if not dm_script.exists():
        print(f"[CRITICAL ERROR] dataset-manager.py not found at: {dm_script}")
        return

    # 2. Configure paths
    base_dir = Path(args.out_dir) 
    uncompressed_root = base_dir / args.genomes_subdir
    base_dir.mkdir(parents=True, exist_ok=True)
    uncompressed_root.mkdir(parents=True, exist_ok=True)

    # 3. Read lists
    selected_items = []
    def collect(file_path):
        if not Path(file_path).exists(): return
        with open(file_path, 'r') as f:
            for line in f:
                path_str = line.strip()
                if not path_str: continue
                acc = normalize_accession(path_str)
                # folder_name is the folder ID used by your pipeline (staphylococcus_..._GCA_...)
                folder_name = Path(path_str).parent.name if "/" in path_str else Path(path_str).name
                if acc:
                    selected_items.append((acc, folder_name))

    collect(args.non_targets)
    collect(args.representatives)

    if not selected_items:
        print("[ERROR] No genomes to process.")
        return

    print(f"[INFO] Processing {len(selected_items)} genomes...")

    final_accessions = []
    index_lines = [] # For kSNP index
    failed = 0

    for i, (acc, folder_name) in enumerate(selected_items, 1):
        zip_dest = base_dir / f"{folder_name}.zip"
        # dataset-manager creates a folder based on zip name or arguments,
        # so we force output to uncompressed/
        
        # A. Download
        if not (zip_dest.exists() and zip_dest.stat().st_size > 5000):
            print(f"[{i}] Downloading {acc}...")
            # IMPORTANT: We include cds, gff3, and seq-report so dataset-manager doesn't fail finding them
            cmd_dl = [datasets_exe, "download", "genome", "accession", acc, 
                      "--include", "genome,seq-report,cds,gff3", 
                      "--filename", str(zip_dest)]
            try:
                subprocess.run(cmd_dl, check=True, capture_output=False)
            except subprocess.CalledProcessError:
                print(f"    [!] NCBI download failed.")
                failed += 1
                continue
        else:
             print(f"[{i}] {acc} already downloaded.")

        # B. Run dataset-manager.py
        # Equivalent command in ksnp-preparation.sh:
        # python3 dataset-manager.py build-dataset --output <dir> <zip>
        print(f"    -> Running dataset-manager...")
        
        cmd_dm = [
            sys.executable, str(dm_script), 
            "build-dataset", 
            "--output", str(uncompressed_root),
            str(zip_dest)
        ]

        try:
            # dataset-manager prints the report line to stdout (Genus Species Strain...)
            # We need to capture it to build the index
            result = subprocess.run(cmd_dm, check=True, capture_output=True, text=True)
            
            # dataset-manager output is a tab-separated line.
            # ksnp-preparation.sh adds the filename at the beginning.
            # DM Output:  Genus \t Species \t Strain \t TaxID \t Path
            # Final Output: FolderName \t Genus \t Species ...
            
            dm_output = result.stdout.strip()
            
            if dm_output:
                # Build kSNP index line
                full_index_line = f"{folder_name}\t{dm_output}"
                index_lines.append(full_index_line)
                final_accessions.append(acc)
                print("    [OK] Processed successfully.")
            else:
                print("    [!] dataset-manager returned no data.")
                failed += 1

        except subprocess.CalledProcessError as e:
            print(f"    [!] Error in dataset-manager: {e.stderr}")
            failed += 1
        except Exception as e:
            print(f"    [!] Unexpected error: {e}")
            failed += 1

    # --- SAVE FILES ---

    # 1. selected_accessions.txt
    print(f"\n[INFO] Saving {args.acc_file}")
    with open(args.acc_file, 'w') as f:
        for acc in sorted(list(set(final_accessions))):
            f.write(f"{acc}\n")

    # 2. ksnp_files_index.tsv (Optional)
    ksnp_index_path = Path(args.acc_file).parent / "ksnp_files_index.tsv"
    print(f"[INFO] Saving kSNP index to {ksnp_index_path}")
    with open(ksnp_index_path, 'w') as f:
        for line in index_lines:
            f.write(f"{line}\n")

    print(f"\n[DONE] Completed: {len(final_accessions)} | Failed: {failed}")

if __name__ == "__main__":
    main()