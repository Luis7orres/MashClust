#!/usr/bin/env python3
"""
Step 0: Download genomes from NCBI (FNA only - With API Key Support)
"""
import subprocess
import sys
import re
import argparse
import json
import gzip
import os
import time
from pathlib import Path
from zipfile import ZipFile

def download_with_adaptive_retry(cmd, batch_num, max_retries=5, api_key=None):
    """
    Download with adaptive retries.
    Uses NCBI API key if provided for higher rate limits.
    """
    base_wait = 60
    
    # Add API key to environment if provided
    env = os.environ.copy()
    if api_key:
        env['NCBI_API_KEY'] = api_key
    
    for attempt in range(max_retries):
        try:
            result = subprocess.run(
                cmd, 
                check=True, 
                capture_output=True, 
                text=True,
                env=env  # Pass environment with API key
            )
            return True, None
        except subprocess.CalledProcessError as e:
            error_msg = e.stderr if e.stderr else str(e)
            
            # Detect error type
            if "giving up after" in error_msg or "gateway" in error_msg.lower():
                error_type = "API_LIMIT"
                wait_multiplier = 2
            elif "INTERNAL_ERROR" in error_msg or "stream error" in error_msg:
                error_type = "STREAM_ERROR"
                wait_multiplier = 1.5
            else:
                error_type = "UNKNOWN"
                wait_multiplier = 1.5
            
            if attempt < max_retries - 1:
                wait_time = int(base_wait * (wait_multiplier ** attempt) + (attempt * 10))
                print(f"[RETRY] Batch {batch_num}, attempt {attempt + 1}/{max_retries} failed ({error_type})", 
                      file=sys.stderr)
                print(f"[RETRY] Waiting {wait_time}s before retry...", file=sys.stderr)
                time.sleep(wait_time)
            else:
                print(f"[FAILED] Batch {batch_num} exhausted all {max_retries} retries", file=sys.stderr)
                return False, error_type
    
    return False, "MAX_RETRIES"

def save_failed_accessions(failed_acc, output_dir, batch_num):
    """Save failed accessions for potential manual retry"""
    failed_file = output_dir / f"failed_batch_{batch_num}.txt"
    failed_file.write_text('\n'.join(failed_acc))
    return failed_file

def main():
    parser = argparse.ArgumentParser(
        description="Step 0: Download genomes in batches with robust error handling"
    )
    parser.add_argument("accessions_file", help="Path to accessions list")
    parser.add_argument("-o", "--output-dir", required=True, help="Output directory")
    parser.add_argument("--keep-zip", action="store_true", help="Keep downloaded zips")
    parser.add_argument("--batch-size", type=int, default=100, 
                       help="Genomes per batch (default: 100 - conservative)")
    parser.add_argument("--delay", type=int, default=95, 
                       help="Delay between successful batches in seconds (default: 90)")
    parser.add_argument("--max-retries", type=int, default=5, 
                       help="Max retries per batch (default: 5)")
    parser.add_argument("--min-success-rate", type=float, default=0.95,
                       help="Minimum success rate to continue (default: 0.95)")
    parser.add_argument("--api-key", type=str, default=None,
                       help="NCBI API key for higher rate limits (RECOMMENDED)")
    parser.add_argument("--start-from-batch", type=int, default=1,
                       help="Resume from specific batch number (default: 1)")
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    temp_dir = out_dir / "temp_downloads"
    genomes_dir = out_dir / "genomes"
    failed_dir = out_dir / "failed_batches"
    
    # Create directories
    temp_dir.mkdir(parents=True, exist_ok=True)
    genomes_dir.mkdir(parents=True, exist_ok=True)
    failed_dir.mkdir(parents=True, exist_ok=True)

    # Load and validate accessions
    accessions = []
    with open(args.accessions_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                if re.match(r'^GC[FA]_\d{9}\.\d+$', line):
                    accessions.append(line)
                else:
                    print(f"[WARNING] Invalid accession format skipped: {line}", file=sys.stderr)

    if not accessions:
        print("[ERROR] No valid accessions found in input file", file=sys.stderr)
        sys.exit(1)

    # Count already downloaded genomes
    existing_genomes = len(list(genomes_dir.rglob("GENOME_*.fna")))
    
    total_genomes = len(accessions)
    total_batches = (total_genomes + args.batch_size - 1) // args.batch_size
    
    print(f"[INFO] ========== MashClust Download Step ==========", file=sys.stderr)
    print(f"[INFO] Total genomes to download: {total_genomes}", file=sys.stderr)
    print(f"[INFO] Already downloaded: {existing_genomes}", file=sys.stderr)
    print(f"[INFO] Batch size: {args.batch_size}", file=sys.stderr)
    print(f"[INFO] Total batches: {total_batches}", file=sys.stderr)
    print(f"[INFO] Starting from batch: {args.start_from_batch}", file=sys.stderr)
    print(f"[INFO] Inter-batch delay: {args.delay}s", file=sys.stderr)
    print(f"[INFO] Max retries per batch: {args.max_retries}", file=sys.stderr)
    print(f"[INFO] Minimum success rate: {args.min_success_rate * 100}%", file=sys.stderr)
    print(f"[INFO] API Key: {'YES (higher limits)' if args.api_key else 'NO (default limits)'}", file=sys.stderr)
    
    if not args.api_key:
        print(f"[WARNING] No API key provided. Consider getting one from:", file=sys.stderr)
        print(f"[WARNING] https://www.ncbi.nlm.nih.gov/account/settings/", file=sys.stderr)
        print(f"[WARNING] Use --api-key YOUR_KEY for 10x higher rate limits", file=sys.stderr)
    
    print(f"[INFO] ============================================", file=sys.stderr)

    extracted_count = existing_genomes
    failed_batches = []
    start_time = time.time()

    # Skip to start batch if resuming
    start_index = (args.start_from_batch - 1) * args.batch_size

    # Process batches
    for i in range(start_index, total_genomes, args.batch_size):
        batch = accessions[i:i + args.batch_size]
        current_batch_num = (i // args.batch_size) + 1
        batch_start_time = time.time()

        print(f"\n[BATCH {current_batch_num}/{total_batches}] Processing {len(batch)} accessions...", 
              file=sys.stderr)

        temp_list_path = temp_dir / f"acc_list_batch_{current_batch_num}.txt"
        temp_list_path.write_text('\n'.join(batch))
        zip_file = temp_dir / f"genomes_batch_{current_batch_num}.zip"

        # Build download command
        cmd = [
            "datasets", "download", "genome", "accession",
            "--inputfile", str(temp_list_path),
            "--filename", str(zip_file)
        ]

        # Execute with adaptive retry
        success, error_type = download_with_adaptive_retry(
            cmd, current_batch_num, args.max_retries, args.api_key
        )
        
        if not success:
            print(f"[ERROR] Batch {current_batch_num} failed after retries (Type: {error_type})", 
                  file=sys.stderr)
            failed_batches.append({
                'batch_num': current_batch_num,
                'error_type': error_type,
                'accessions': batch,
                'count': len(batch)
            })
            
            # Save failed accessions
            failed_file = save_failed_accessions(batch, failed_dir, current_batch_num)
            print(f"[INFO] Failed accessions saved: {failed_file}", file=sys.stderr)
            
            # Clean temp files
            if temp_list_path.exists():
                temp_list_path.unlink()
            if zip_file.exists():
                zip_file.unlink()
            
            # Check if we should continue based on success rate
            total_processed = i + len(batch)
            current_success_rate = extracted_count / total_processed if total_processed > 0 else 0
            
            if current_success_rate < args.min_success_rate and current_batch_num > 10:
                print(f"[FATAL] Success rate ({current_success_rate:.2%}) below minimum ({args.min_success_rate:.2%})", 
                      file=sys.stderr)
                print(f"[FATAL] Too many consecutive failures - stopping pipeline", file=sys.stderr)
                print(f"[SUGGESTION] Wait 2 hours and resume with: --start-from-batch {current_batch_num}", file=sys.stderr)
                sys.exit(1)
            
            continue

        # Extract and organize current batch
        batch_extracted = 0
        try:
            if not zip_file.exists():
                raise FileNotFoundError(f"Downloaded zip file not found: {zip_file}")
                
            with ZipFile(zip_file, 'r') as z:
                # 1. Read batch metadata
                genome_meta = {}
                for filename in z.namelist():
                    if filename.endswith('assembly_data_report.jsonl'):
                        with z.open(filename) as f:
                            for line in f:
                                data = json.loads(line.decode('utf-8'))
                                acc = data.get('accession')
                                if acc:
                                    genome_meta[acc] = {
                                        'organism': data.get('organism', {}).get('organismName', 'unknown'),
                                        'tax_id': data.get('organism', {}).get('taxId', 'unknown')
                                    }

                # 2. Extract .fna files from batch
                for filename in z.namelist():
                    if filename.endswith(('.fna', '.fna.gz')) and '/data/GC' in filename:
                        match = re.search(r'(GC[FA]_\d{9}\.\d+)', filename)
                        if match:
                            acc = match.group(1)
                            if acc in batch:
                                meta = genome_meta.get(acc, {})
                                org_name = re.sub(
                                    r'[^a-zA-Z0-9-_]', '', 
                                    meta.get('organism', 'unknown').replace(' ', '_').lower()
                                )
                                tax_id = meta.get('tax_id', 'unknown')
                                acc_clean = acc.replace('.', '_')

                                dest_folder = genomes_dir / f"{org_name}_{tax_id}_{acc_clean}"
                                dest_folder.mkdir(exist_ok=True)
                                dest_file = dest_folder / f"GENOME_{acc_clean}.fna"

                                with z.open(filename) as source, open(dest_file, 'wb') as target:
                                    if filename.endswith('.gz'):
                                        target.write(gzip.decompress(source.read()))
                                    else:
                                        target.write(source.read())
                                batch_extracted += 1
                                extracted_count += 1

            batch_time = time.time() - batch_start_time
            print(f"[SUCCESS] Batch {current_batch_num}: {batch_extracted} genomes in {batch_time:.1f}s", 
                  file=sys.stderr)
            print(f"[PROGRESS] Total: {extracted_count}/{total_genomes} ({extracted_count/total_genomes*100:.1f}%)", 
                  file=sys.stderr)

        except Exception as e:
            print(f"[ERROR] Extraction failed for batch {current_batch_num}: {e}", file=sys.stderr)
            failed_batches.append({
                'batch_num': current_batch_num,
                'error_type': 'EXTRACTION_ERROR',
                'accessions': batch,
                'count': len(batch)
            })
            save_failed_accessions(batch, failed_dir, current_batch_num)
            
        finally:
            # Clean temporary files
            if not args.keep_zip and zip_file.exists():
                zip_file.unlink()
            if temp_list_path.exists():
                temp_list_path.unlink()

        # Wait before next batch
        if current_batch_num < total_batches:
            print(f"[WAIT] Pausing {args.delay}s before next batch...", file=sys.stderr)
            time.sleep(args.delay)

    # Final statistics and genome list generation
    elapsed_time = time.time() - start_time
    elapsed_hours = elapsed_time / 3600
    success_rate = (extracted_count / total_genomes) * 100 if total_genomes > 0 else 0
    failed_count = sum(fb['count'] for fb in failed_batches)
    
    print(f"\n{'='*70}", file=sys.stderr)
    print(f"[SUMMARY] Download Process Completed", file=sys.stderr)
    print(f"{'='*70}", file=sys.stderr)
    print(f"  Total genomes requested:    {total_genomes}", file=sys.stderr)
    print(f"  Successfully extracted:     {extracted_count}", file=sys.stderr)
    print(f"  Failed:                     {failed_count}", file=sys.stderr)
    print(f"  Success rate:               {success_rate:.2f}%", file=sys.stderr)
    print(f"  Total time:                 {elapsed_hours:.2f} hours", file=sys.stderr)
    
    if failed_batches:
        print(f"\n[FAILED BATCHES]", file=sys.stderr)
        for fb in failed_batches:
            print(f"  Batch {fb['batch_num']:3d}: {fb['error_type']:20s} ({fb['count']} genomes)", 
                  file=sys.stderr)
        print(f"\n[INFO] Failed accessions directory: {failed_dir}", file=sys.stderr)
    
    print(f"{'='*70}", file=sys.stderr)

    # Generate genome list for Mash
    print(f"\n[INFO] Generating genome list for Mash...", file=sys.stderr)
    genome_list_file = out_dir / "genome_list.txt"
    genome_count = 0
    
    with open(genome_list_file, "w") as f:
        for p in sorted(genomes_dir.rglob("GENOME_*.fna")):
            f.write(f"{p}\n")
            genome_count += 1
    
    print(f"[INFO] Genome list created: {genome_list_file}", file=sys.stderr)
    print(f"[INFO] Total genomes in list: {genome_count}", file=sys.stderr)

    if genome_count == 0:
        print(f"[FATAL] No genomes were successfully downloaded", file=sys.stderr)
        sys.exit(1)
    
    if success_rate < (args.min_success_rate * 100):
        print(f"[WARNING] Success rate ({success_rate:.2f}%) below target", file=sys.stderr)
    
    if failed_batches:
        print(f"\n[WARNING] Pipeline completed with {len(failed_batches)} failed batches", file=sys.stderr)
    else:
        print(f"\n[SUCCESS] All batches completed successfully!", file=sys.stderr)
    
    sys.exit(0)

if __name__ == "__main__":
    main()