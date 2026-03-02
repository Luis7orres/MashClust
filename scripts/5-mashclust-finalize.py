#!/usr/bin/env python
"""
Script 5: Finalize Results and Integrate with dataset-manager.py
Restored with --genomes-subdir support while maintaining new features.
"""
import argparse
import shutil
import re
import os
import sys
import subprocess
import time
import zipfile
from pathlib import Path

def normalize_accession(path_str):
    """
    Extracts GCA_000000000.1 from strings like paths or folder names.
    Returns None if no valid accession is found.
    """
    pattern = r'(GC[FA])_(\d{9})[_.](\d+)'
    match = re.search(pattern, path_str)
    if match:
        return f"{match.group(1)}_{match.group(2)}.{match.group(3)}"
    return None

def check_zip_integrity(zip_path):
    """
    Validates if the ZIP file exists and passes the CRC-32 consistency check.
    Returns True if the file is healthy, False otherwise.
    """
    if not zip_path.exists():
        return False
    try:
        with zipfile.ZipFile(zip_path) as zf:
            # testzip() returns the name of the first corrupt file or None if all OK
            return zf.testzip() is None
    except Exception:
        return False

def main():
    parser = argparse.ArgumentParser(description="Step 5: Finalize Results")
    parser.add_argument("--non-targets", required=True, help="List of non-target genome paths")
    parser.add_argument("--representatives", required=True, help="List of representative genome paths")
    parser.add_argument("--out-dir", required=True, help="Directory to save the results")
    parser.add_argument("--acc-file", required=True, help="Path for the selected_accessions.txt output")
    parser.add_argument("--genomes-subdir", default="uncompressed", help="Subdirectory for uncompressed genomes")
    parser.add_argument("--dataset-manager", required=True, help="Absolute path to dataset-manager.py")
    parser.add_argument("--api-key", default=None, help="NCBI API Key for higher rate limits")
    args = parser.parse_args()

    # Tool and path verification
    datasets_exe = shutil.which("datasets")
    if not datasets_exe:
        print("[CRITICAL ERROR] 'datasets' binary not found in PATH.")
        sys.exit(1)

    dm_script = Path(args.dataset_manager)
    if not dm_script.exists():
        print(f"[CRITICAL ERROR] dataset-manager.py not found at: {dm_script}")
        sys.exit(1)

    # Output directory setup
    base_dir = Path(args.out_dir) 
    uncompressed_root = base_dir / args.genomes_subdir
    base_dir.mkdir(parents=True, exist_ok=True)
    uncompressed_root.mkdir(parents=True, exist_ok=True)

    # Collect genomes to process
    selected_items = []
    def collect(file_path):
        if not Path(file_path).exists(): return
        with open(file_path, 'r') as f:
            for line in f:
                path_str = line.strip()
                if not path_str: continue
                acc = normalize_accession(path_str)
                folder_name = Path(path_str).parent.name if "/" in path_str else Path(path_str).name
                if acc:
                    selected_items.append((acc, folder_name))

    collect(args.non_targets)
    collect(args.representatives)

    if not selected_items:
        print("[ERROR] No accessions found to process.")
        return

    print(f"[INFO] Processing {len(selected_items)} genomes...")
    final_accessions = []
    index_lines = [] 
    failed = 0

    # Configure environment with API Key
    env = os.environ.copy()
    if args.api_key:
        env['NCBI_API_KEY'] = args.api_key

    for i, (acc, folder_name) in enumerate(selected_items, 1):
        zip_dest = base_dir / f"{folder_name}.zip"
        success_dl = False
        
        # DOWNLOAD LOOP WITH 5 RETRIES AND INTEGRITY CHECK
        for attempt in range(1, 6): # 5 attempts
            if check_zip_integrity(zip_dest):
                success_dl = True
                if attempt == 1: print(f"[{i}] {acc} already exists and is valid.")
                break
            else:
                if zip_dest.exists():
                    print(f"    [!] Corrupt ZIP detected (CRC-32 failure) for {acc}. Deleting and retrying...")
                    zip_dest.unlink()

                print(f"[{i}] Downloading {acc} (Attempt {attempt}/5)...")
                cmd_dl = [
                    datasets_exe, "download", "genome", "accession", acc, 
                    "--include", "genome,seq-report,cds,gff3", 
                    "--filename", str(zip_dest)
                ]
                try:
                    subprocess.run(cmd_dl, check=True, capture_output=True, env=env)
                    if check_zip_integrity(zip_dest):
                        success_dl = True
                        break
                except subprocess.CalledProcessError:
                    pass
            
            if attempt < 5:
                time.sleep(2)

        if not success_dl:
            print(f"    [FATAL] NCBI download failed permanently for {acc} after 5 attempts.")
            failed += 1
            continue

        # RUN DATASET-MANAGER.PY
        print(f"    -> Running dataset-manager for extraction...")
        cmd_dm = [
            sys.executable, str(dm_script), 
            "build-dataset", 
            "--output", str(uncompressed_root),
            str(zip_dest)
        ]

        try:
            result = subprocess.run(cmd_dm, check=True, capture_output=True, text=True)
            dm_output = result.stdout.strip()
            
            if dm_output:
                index_lines.append(f"{folder_name}\t{dm_output}")
                final_accessions.append(acc)
                print("    [OK] Extracted successfully.")
            else:
                print("    [!] dataset-manager returned no output.")
                failed += 1

        except subprocess.CalledProcessError as e:
            print(f"    [!] dataset-manager failed: {e.stderr}")
            failed += 1

    # Save final lists
    print(f"\n[INFO] Saving accession list to {args.acc_file}")
    with open(args.acc_file, 'w') as f:
        for acc in sorted(list(set(final_accessions))):
            f.write(f"{acc}\n")

    ksnp_index_path = Path(args.acc_file).parent / "ksnp_files_index.tsv"
    print(f"[INFO] Saving kSNP index to {ksnp_index_path}")
    with open(ksnp_index_path, 'w') as f:
        for line in index_lines:
            f.write(f"{line}\n")

    print(f"\n[DONE] Completed: {len(final_accessions)} | Failed: {failed}")

if __name__ == "__main__":
    main()