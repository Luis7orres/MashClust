#!/usr/bin/env python3
"""
1-checkm-filter.py - MARGE Quality Filter
Logic:
1. Flag 'suppressed' assemblies.
2. Flag genomes below completeness threshold.
3. Flag genomes with missing CheckM data.
All flagged genomes are added to the blacklist.
If --clean-dir is provided, flagged genome folders are physically deleted.
"""

import json
import argparse
import sys
import shutil
from pathlib import Path


# ── helpers ───────────────────────────────────────────────────────────────────

def accession_with_underscore(accession: str) -> str:
    """Standardize accession format for folder matching (GCA_001.1 -> GCA_001_1)."""
    return accession.replace(".", "_")


def parse_report(record: dict) -> dict:
    """
    Extract metadata from NCBI JSONL records.
    Tries both snake_case and camelCase field names for robustness.
    """
    acc = (record.get("current_accession")
           or record.get("currentAccession")
           or record.get("accession", "UNKNOWN"))

    organism = record.get("organism", {})
    org = organism.get("organism_name") or organism.get("organismName", "N/A")

    ainfo = record.get("assembly_info") or record.get("assemblyInfo", {})
    level  = ainfo.get("assembly_level")  or ainfo.get("assemblyLevel",  "N/A")
    status = ainfo.get("assembly_status") or ainfo.get("assemblyStatus", "N/A")
    tech   = ainfo.get("sequencing_tech") or ainfo.get("sequencingTech", "N/A")
    notes  = ainfo.get("genome_notes")    or ainfo.get("genomeNotes",    [])

    checkm = record.get("checkm_info") or record.get("checkmInfo", {})
    comp   = checkm.get("completeness")
    cont   = checkm.get("contamination")
    mset   = (checkm.get("checkm_marker_set")
              or checkm.get("markerSet")
              or checkm.get("marker_set", "N/A"))

    return {
        "accession":       acc,
        "organism":        org,
        "assembly_level":  level,
        "status":          status,
        "sequencing_tech": tech,
        "genome_notes":    notes if isinstance(notes, list) else [],
        "completeness":    comp,   # None if missing
        "contamination":   cont,   # None if missing
        "marker_set":      mset,
    }


def load_jsonl_file(path: Path) -> list[dict]:
    """
    Accepts both:
      - true JSONL  (one JSON object per line)
      - a single JSON object / array spanning multiple lines
    Returns a flat list of raw record dicts.
    """
    records = []
    raw = path.read_text(encoding="utf-8").strip()

    # Try line-by-line JSONL first
    lines = [l.strip() for l in raw.splitlines() if l.strip()]
    parsed_lines = []
    for line in lines:
        try:
            parsed_lines.append(json.loads(line))
        except json.JSONDecodeError:
            parsed_lines = []
            break

    if parsed_lines:
        objects = parsed_lines
    else:
        try:
            objects = [json.loads(raw)]
        except json.JSONDecodeError as e:
            print(f"  [WARN] Cannot parse {path.name}: {e}", file=sys.stderr)
            return []

    for obj in objects:
        if isinstance(obj, dict) and "reports" in obj:
            records.extend(obj["reports"])
        elif isinstance(obj, list):
            records.extend(obj)
        elif isinstance(obj, dict):
            records.append(obj)

    return records


# ── report builder ────────────────────────────────────────────────────────────

def build_report(passing, failing, suppressed, no_checkm, threshold, folder) -> str:
    SEP  = "=" * 72
    DASH = "─" * 72
    total = len(passing) + len(failing) + len(suppressed) + len(no_checkm)
    lines = []

    lines += [
        SEP,
        "  MARGE GENOME COMPLETENESS FILTER REPORT",
        SEP,
        f"  Folder          : {Path(folder).resolve()}",
        f"  Threshold       : completeness >= {threshold}%",
        f"  Total records   : {total}",
        f"  Passing         : {len(passing)}",
        f"  Failing (<{threshold}%) : {len(failing)}",
        f"  Suppressed      : {len(suppressed)}",
        f"  No CheckM data  : {len(no_checkm)}",
        f"  Blacklisted     : {len(failing) + len(suppressed) + len(no_checkm)}",
        SEP,
    ]

    if suppressed:
        lines += ["", DASH, "  SUPPRESSED ASSEMBLIES", DASH]
        for i, g in enumerate(suppressed, 1):
            lines.append(
                f"\n  [{i}] {g['accession']}"
                f"\n      Organism        : {g['organism']}"
                f"\n      Assembly level  : {g['assembly_level']}"
                f"\n      Assembly status : {g['status']}"
            )

    if failing:
        lines += ["", DASH, f"  FAILING GENOMES  (completeness < {threshold}%)", DASH]
        for i, g in enumerate(failing, 1):
            notes_str = "; ".join(g["genome_notes"]) if g["genome_notes"] else "none"
            comp_str  = f"{g['completeness']:.2f}%" if g["completeness"] is not None else "N/A"
            cont_str  = f"{g['contamination']:.2f}%" if g["contamination"] is not None else "N/A"
            lines.append(
                f"\n  [{i}] {g['accession']}"
                f"\n      Organism        : {g['organism']}"
                f"\n      Assembly level  : {g['assembly_level']}"
                f"\n      Assembly status : {g['status']}"
                f"\n      Sequencing tech : {g['sequencing_tech']}"
                f"\n      Completeness    : {comp_str}"
                f"\n      Contamination   : {cont_str}"
                f"\n      CheckM marker   : {g['marker_set']}"
                f"\n      Genome notes    : {notes_str}"
            )

    if no_checkm:
        lines += ["", DASH, "  RECORDS WITHOUT CheckM DATA  (blacklisted)", DASH]
        for g in no_checkm:
            lines.append(f"  {g['accession']}  |  {g['organism']}")

    if passing:
        lines += ["", DASH, "  PASSING GENOMES", DASH]
        for g in passing:
            comp_str = f"{g['completeness']:.2f}%" if g["completeness"] is not None else "N/A"
            lines.append(f"  {g['accession']}  |  {comp_str}  |  {g['organism']}")

    lines.append(f"\n{SEP}\n")
    return "\n".join(lines)


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Filter out poor quality or suppressed genomes."
    )
    parser.add_argument("json_folder",  help="Folder with .jsonl / .json reports")
    parser.add_argument("--threshold",  type=float, default=90.0,
                        help="Min completeness %% (default: 90.0)")
    parser.add_argument("--output",     required=True,
                        help="Report file path (pipeline expects quality_report.txt)")
    parser.add_argument("--clean-dir",
                        help="Optional: genome folder — flagged subfolders will be deleted")
    # --- NUEVO: Argumento de la whitelist ---
    parser.add_argument("--whitelist", nargs='*', default=[],
                        help="Accessions to always keep regardless of quality")
    args = parser.parse_args()

    json_dir = Path(args.json_folder)
    if not json_dir.exists():
        print(f"[ERROR] Reports directory not found: {json_dir}", file=sys.stderr)
        sys.exit(1)

    files = sorted(json_dir.glob("*.jsonl")) + sorted(json_dir.glob("*.json"))
    if not files:
        print(f"[ERROR] No .json/.jsonl files found in '{json_dir}'", file=sys.stderr)
        sys.exit(1)

    passing, failing, suppressed, no_checkm = [], [], [], []
    
    # --- NUEVO: Convertimos la lista a un set para búsquedas rápidas ---
    whitelist = set(args.whitelist) if args.whitelist else set()

    # 1. Parse all report files
    for jfile in files:
        for rec in load_jsonl_file(jfile):
            g = parse_report(rec)

            # --- NUEVO: Si está en la whitelist, lo pasamos directamente y saltamos el resto de comprobaciones ---
            if g["accession"] in whitelist:
                passing.append(g)
                continue

            if g["status"].lower() == "suppressed":
                suppressed.append(g)
            elif g["completeness"] is None:
                no_checkm.append(g)
            elif g["completeness"] < args.threshold:
                failing.append(g)
            else:
                passing.append(g)

    # 2. Write rich report
    report_path = Path(args.output)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_text = build_report(passing, failing, suppressed, no_checkm,
                               args.threshold, args.json_folder)
    report_path.write_text(report_text, encoding="utf-8")

    # 3. Write accession blacklist + optional physical deletion
    blacklist = failing + suppressed + no_checkm
    acc_path  = report_path.with_name(report_path.stem + "_accessions.txt")

    with open(acc_path, "w") as f_acc:
        for g in blacklist:
            acc_id = accession_with_underscore(g["accession"])
            f_acc.write(f"{acc_id}\n")

            if args.clean_dir:
                clean_path = Path(args.clean_dir)
                for folder in clean_path.glob(f"*{acc_id}*"):
                    if folder.is_dir():
                        shutil.rmtree(folder)
                        print(f"[CLEANUP] Deleted: {folder.name}")

    total = len(passing) + len(blacklist)
    print(f"[INFO] Processed {total} genomes.")
    print(f"[INFO] Passing: {len(passing)} | Failing: {len(failing)} | "
          f"Suppressed: {len(suppressed)} | No CheckM: {len(no_checkm)}")
    print(f"[INFO] Report    : {report_path}")
    print(f"[INFO] Blacklist : {acc_path}  ({len(blacklist)} entries)")


if __name__ == "__main__":
    main()