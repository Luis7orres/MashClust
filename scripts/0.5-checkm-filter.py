#!/usr/bin/env python3
"""
0.5-checkm-filter.py - MARGE Quality Filter
Logic:
  Each of the three filter stages is independently togglable:
    --filter-suppressed   : blacklist assemblies with status 'suppressed'
    --filter-no-checkm    : blacklist records missing CheckM metadata
    --filter-completeness : blacklist genomes below --threshold
  Any combination is valid. All flagged genomes go to the blacklist.
  If --clean-dir is provided, flagged genome folders are physically deleted.
"""

import json
import argparse
import sys
import shutil
from pathlib import Path


# ── helpers ───────────────────────────────────────────────────────────────────

def accession_with_underscore(accession: str) -> str:
    return accession.replace(".", "_")


def parse_report(record: dict) -> dict:
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
        "completeness":    comp,
        "contamination":   cont,
        "marker_set":      mset,
    }


def load_jsonl_file(path: Path) -> list[dict]:
    records = []
    raw = path.read_text(encoding="utf-8").strip()

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

def build_report(passing, failing, suppressed, no_checkm,
                 threshold, folder, filters_active: dict,
                 blacklist_count: int) -> str:
    SEP  = "=" * 72
    DASH = "─" * 72
    total = len(passing) + len(failing) + len(suppressed) + len(no_checkm)

    active_labels = []
    if filters_active["suppressed"]:
        active_labels.append("suppressed")
    if filters_active["no_checkm"]:
        active_labels.append("missing CheckM")
    if filters_active["completeness"]:
        active_labels.append(f"completeness < {threshold}%")
    active_str = ", ".join(active_labels) if active_labels else "none"

    lines = [
        SEP,
        "  MARGE GENOME COMPLETENESS FILTER REPORT",
        SEP,
        f"  Folder          : {Path(folder).resolve()}",
        f"  Active filters  : {active_str}",
        f"  Threshold       : completeness >= {threshold}%",
        f"  Total records   : {total}",
        f"  Passing         : {len(passing)}",
        f"  Failing (<{threshold}%) : {len(failing) if filters_active['completeness'] else 0}",
        f"  Suppressed      : {len(suppressed) if filters_active['suppressed'] else 0}",
        f"  No CheckM data  : {len(no_checkm) if filters_active['no_checkm'] else 0}",
        f"  Blacklisted     : {blacklist_count}",
        SEP,
    ]

    if suppressed:
        label = "SUPPRESSED ASSEMBLIES (BLACKLISTED)" if filters_active["suppressed"] \
                else "SUPPRESSED ASSEMBLIES (filter OFF — kept)"
        lines += ["", DASH, f"  {label}", DASH]
        for i, g in enumerate(suppressed, 1):
            lines.append(
                f"\n  [{i}] {g['accession']}"
                f"\n      Organism        : {g['organism']}"
                f"\n      Assembly level  : {g['assembly_level']}"
                f"\n      Assembly status : {g['status']}"
            )

    if failing:
        label = f"FAILING GENOMES  (completeness < {threshold}%)" if filters_active["completeness"] \
                else f"FAILING GENOMES  (completeness < {threshold}%) (filter OFF — kept)"
        lines += ["", DASH, f"  {label}", DASH]
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
        label = "RECORDS WITHOUT CheckM DATA (BLACKLISTED)" if filters_active["no_checkm"] \
                else "RECORDS WITHOUT CheckM DATA (filter OFF — kept)"
        lines += ["", DASH, f"  {label}", DASH]
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
        description="Filter poor-quality or suppressed genomes. Each stage is independently togglable."
    )
    parser.add_argument("json_folder",  help="Folder with .jsonl / .json reports")
    parser.add_argument("--threshold",  type=float, default=90.0,
                        help="Min completeness %% (default: 90.0)")
    parser.add_argument("--output",     required=True,
                        help="Report file path (pipeline expects quality_report.txt)")
    parser.add_argument("--clean-dir",
                        help="Optional: genome folder — flagged subfolders will be deleted")
    parser.add_argument("--whitelist",  nargs="*", default=[],
                        help="Accessions to always keep regardless of quality")

    # ── per-stage filter toggles ──────────────────────────────────────────────
    parser.add_argument("--filter-suppressed",   dest="filter_suppressed",
                        action=argparse.BooleanOptionalAction, default=True,
                        help="Blacklist assemblies with status 'suppressed' (default: on)")
    parser.add_argument("--filter-no-checkm",    dest="filter_no_checkm",
                        action=argparse.BooleanOptionalAction, default=True,
                        help="Blacklist records missing CheckM metadata (default: on)")
    parser.add_argument("--filter-completeness", dest="filter_completeness",
                        action=argparse.BooleanOptionalAction, default=True,
                        help="Blacklist genomes below --threshold (default: on)")

    args = parser.parse_args()

    filters_active = {
        "suppressed":   args.filter_suppressed,
        "no_checkm":    args.filter_no_checkm,
        "completeness": args.filter_completeness,
    }

    json_dir = Path(args.json_folder)
    if not json_dir.exists():
        print(f"[ERROR] Reports directory not found: {json_dir}", file=sys.stderr)
        sys.exit(1)

    files = sorted(json_dir.glob("*.jsonl")) + sorted(json_dir.glob("*.json"))
    if not files:
        print(f"[ERROR] No .json/.jsonl files found in '{json_dir}'", file=sys.stderr)
        sys.exit(1)

    passing, failing, suppressed, no_checkm = [], [], [], []
    whitelist = set(args.whitelist)

    for jfile in files:
        for rec in load_jsonl_file(jfile):
            g = parse_report(rec)

            # Whitelist always wins
            if g["accession"] in whitelist:
                passing.append(g)
                continue

            is_suppressed = g["status"].lower() == "suppressed"
            is_no_checkm  = g["completeness"] is None
            is_failing    = (not is_suppressed
                             and not is_no_checkm
                             and g["completeness"] < args.threshold)

            if is_suppressed:
                suppressed.append(g)
            elif is_no_checkm:
                no_checkm.append(g)
            elif is_failing:
                failing.append(g)
            else:
                passing.append(g)

    # Blacklist sólo incluye los grupos con filtro activo
    blacklist = []
    if filters_active["suppressed"]:
        blacklist.extend(suppressed)
    else:
        passing.extend(suppressed)

    if filters_active["no_checkm"]:
        blacklist.extend(no_checkm)
    else:
        passing.extend(no_checkm)

    if filters_active["completeness"]:
        blacklist.extend(failing)
    else:
        passing.extend(failing)

    # Write rich report
    report_path = Path(args.output)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_text = build_report(
        passing,
        failing,
        suppressed,
        no_checkm,
        args.threshold,
        args.json_folder,
        filters_active,
        len(blacklist),   # ← contador real
    )
    report_path.write_text(report_text, encoding="utf-8")

    # Write blacklist + optional physical deletion
    acc_path = report_path.with_name(report_path.stem + "_accessions.txt")
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
    print(f"[INFO] Active filters : suppressed={filters_active['suppressed']}  "
          f"no_checkm={filters_active['no_checkm']}  "
          f"completeness={filters_active['completeness']}")
    print(f"[INFO] Passing: {len(passing)} | Failing: {len(failing) if filters_active['completeness'] else 0} | "
          f"Suppressed: {len(suppressed) if filters_active['suppressed'] else 0} | "
          f"No CheckM: {len(no_checkm) if filters_active['no_checkm'] else 0}")
    print(f"[INFO] Report    : {report_path}")
    print(f"[INFO] Blacklist : {acc_path}  ({len(blacklist)} entries)")


if __name__ == "__main__":
    main()