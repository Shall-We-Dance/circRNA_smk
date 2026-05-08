import csv
import math
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path


def parse_float(value):
    try:
        parsed = float(str(value).strip())
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def safe_name(value):
    allowed = set("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789._-")
    cleaned = "".join(ch if ch in allowed else "_" for ch in str(value))
    return cleaned.strip("._") or "event"


def row_is_significant(row, fdr_cutoff, pvalue_cutoff, inc_diff_cutoff):
    fdr = parse_float(row.get("FDR"))
    pvalue = parse_float(row.get("PValue"))
    inc_diff = parse_float(row.get("IncLevelDifference"))

    if fdr is None or inc_diff is None:
        return False
    if fdr > fdr_cutoff:
        return False
    if pvalue is not None and pvalue > pvalue_cutoff:
        return False
    return abs(inc_diff) >= inc_diff_cutoff


def sort_key(row):
    fdr = parse_float(row.get("FDR"))
    pvalue = parse_float(row.get("PValue"))
    inc_diff = parse_float(row.get("IncLevelDifference"))
    return (
        fdr if fdr is not None else float("inf"),
        pvalue if pvalue is not None else float("inf"),
        -(abs(inc_diff) if inc_diff is not None else 0),
    )


event_file = Path(snakemake.input.event_file)
outdir = Path(snakemake.params.outdir)
events_dir = outdir / "filtered_events"
plots_root = outdir / "plots"
manifest_path = Path(snakemake.output.manifest)
done_path = Path(snakemake.output.done)
log_path = Path(snakemake.log[0])
event_type = str(snakemake.wildcards.event_type)

if outdir.exists():
    shutil.rmtree(outdir)
events_dir.mkdir(parents=True, exist_ok=True)
plots_root.mkdir(parents=True, exist_ok=True)
manifest_path.parent.mkdir(parents=True, exist_ok=True)
done_path.parent.mkdir(parents=True, exist_ok=True)
log_path.parent.mkdir(parents=True, exist_ok=True)

fdr_cutoff = float(snakemake.params.fdr_cutoff)
pvalue_cutoff = float(snakemake.params.pvalue_cutoff)
inc_diff_cutoff = float(snakemake.params.inc_diff_cutoff)
max_events = int(snakemake.params.max_events)
fail_on_error = bool(snakemake.params.fail_on_error)

with open(event_file) as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    if reader.fieldnames is None:
        raise ValueError(f"rMATS event file has no header: {event_file}")
    fieldnames = reader.fieldnames
    rows = [
        row
        for row in reader
        if row_is_significant(row, fdr_cutoff, pvalue_cutoff, inc_diff_cutoff)
    ]

rows.sort(key=sort_key)
if max_events > 0:
    rows = rows[:max_events]

base_cmd = [
    sys.executable,
    str(Path(snakemake.input.script).resolve()),
    "--b1",
    str(Path(snakemake.input.b1).resolve()),
    "--b2",
    str(Path(snakemake.input.b2).resolve()),
    "--event-type",
    event_type,
    "--l1",
    str(snakemake.params.label_b1),
    "--l2",
    str(snakemake.params.label_b2),
    "--exon_s",
    str(snakemake.params.exon_s),
    "--intron_s",
    str(snakemake.params.intron_s),
    "--group-info",
    str(Path(snakemake.input.group_info).resolve()),
    "--min-counts",
    str(snakemake.params.min_counts),
    "--font-size",
    str(snakemake.params.font_size),
    "--fig-width",
    str(snakemake.params.fig_width),
]
if float(snakemake.params.fig_height) > 0:
    base_cmd.extend(["--fig-height", str(snakemake.params.fig_height)])
colors = list(snakemake.params.colors)
if colors:
    base_cmd.extend(["--color", ",".join(colors)])
if bool(snakemake.params.keep_event_chr_prefix):
    base_cmd.append("--keep-event-chr-prefix")
if bool(snakemake.params.remove_event_chr_prefix):
    base_cmd.append("--remove-event-chr-prefix")
extra_args = str(snakemake.params.extra_args or "")
if extra_args:
    base_cmd.extend(shlex.split(extra_args))

env = os.environ.copy()
env.setdefault("MPLBACKEND", "Agg")
manifest_rows = []

with open(log_path, "w") as log_handle:
    log_handle.write(
        f"Input event file: {event_file}\n"
        f"Event type: {event_type}\n"
        f"Selected events: {len(rows)}\n"
        f"Thresholds: FDR <= {fdr_cutoff}, PValue <= {pvalue_cutoff}, "
        f"|IncLevelDifference| >= {inc_diff_cutoff}\n\n"
    )

    for index, row in enumerate(rows, start=1):
        event_id = row.get("ID", str(index))
        gene = row.get("geneSymbol") or row.get("GeneID") or "gene"
        event_key = f"{index:05d}_{event_type}_{safe_name(gene)}_{safe_name(event_id)}"
        single_event_file = events_dir / f"{event_key}.txt"
        plot_outdir = plots_root / event_key
        plot_outdir.mkdir(parents=True, exist_ok=True)

        with open(single_event_file, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerow(row)

        cmd = base_cmd + [
            "-e",
            str(single_event_file.resolve()),
            "-o",
            str(plot_outdir.resolve()),
        ]
        log_handle.write(f"[{event_key}] command: {' '.join(shlex.quote(part) for part in cmd)}\n")
        proc = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=env,
        )
        if proc.stdout:
            log_handle.write(proc.stdout)
        if proc.stderr:
            log_handle.write(proc.stderr)

        pdfs = sorted(str(path.resolve()) for path in (plot_outdir / "Sashimi_plot").glob("*.pdf"))
        status = "ok" if proc.returncode == 0 and pdfs else "failed"
        if status == "failed" and fail_on_error:
            raise RuntimeError(
                f"rmats2sashimiplot failed for {event_key}; see {log_path}"
            )
        manifest_rows.append(
            {
                "comparison": snakemake.wildcards.comparison,
                "event_type": event_type,
                "count_type": snakemake.wildcards.count_type,
                "rank": index,
                "event_id": event_id,
                "gene_id": row.get("GeneID", ""),
                "gene_symbol": row.get("geneSymbol", ""),
                "fdr": row.get("FDR", ""),
                "pvalue": row.get("PValue", ""),
                "inc_level_difference": row.get("IncLevelDifference", ""),
                "event_file": str(single_event_file.resolve()),
                "plot_dir": str(plot_outdir.resolve()),
                "pdfs": ",".join(pdfs),
                "status": status,
                "returncode": proc.returncode,
            }
        )
        log_handle.write(f"[{event_key}] status={status} pdfs={len(pdfs)} returncode={proc.returncode}\n\n")

with open(manifest_path, "w", newline="") as handle:
    fieldnames = [
        "comparison",
        "event_type",
        "count_type",
        "rank",
        "event_id",
        "gene_id",
        "gene_symbol",
        "fdr",
        "pvalue",
        "inc_level_difference",
        "event_file",
        "plot_dir",
        "pdfs",
        "status",
        "returncode",
    ]
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(manifest_rows)

done_path.write_text(
    f"events_selected\t{len(rows)}\n"
    f"plots_ok\t{sum(1 for row in manifest_rows if row['status'] == 'ok')}\n"
    f"plots_failed\t{sum(1 for row in manifest_rows if row['status'] == 'failed')}\n"
)
