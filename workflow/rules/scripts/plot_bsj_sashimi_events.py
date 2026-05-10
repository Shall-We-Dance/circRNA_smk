import csv
import math
import os
import re
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


def unique_names(names):
    seen = {}
    out = []
    for name in names:
        base = name or "column"
        if base not in seen:
            seen[base] = 0
            out.append(base)
        else:
            seen[base] += 1
            out.append(f"{base}_{seen[base]}")
    return out


def read_tabular_loose(path):
    lines = [
        line.rstrip("\n")
        for line in Path(path).read_text().splitlines()
        if line.strip()
    ]
    if not lines:
        return [], []

    split_lines = [line.split("\t") for line in lines]
    header = split_lines[0]
    rows = split_lines[1:]
    if not rows:
        return unique_names(header), []

    max_cols = max(len(header), *(len(row) for row in rows))
    if len(header) == max_cols - 1:
        header = ["row_id"] + header
    elif len(header) < max_cols:
        header = header + [f"extra_{idx}" for idx in range(1, max_cols - len(header) + 1)]
    elif len(header) > max_cols:
        header = header[:max_cols]
    header = unique_names(header)

    out_rows = []
    for row in rows:
        padded = list(row)
        padded.extend([""] * (max_cols - len(padded)))
        out_rows.append(dict(zip(header, padded[:max_cols])))
    return header, out_rows


def first_existing(row, candidates):
    for candidate in candidates:
        if candidate in row:
            value = str(row.get(candidate, "")).strip()
            if value:
                return candidate, value
    return None, ""


def looks_like_circ(value):
    return re.match(r"^[^:]+:\d+(?:\||-|\.\.)\d+", str(value).strip()) is not None


def split_circ_gene(value):
    raw = str(value).strip()
    circ = raw.split(";", 1)[0]
    gene = raw.split(";", 1)[1] if ";" in raw else ""
    return circ, gene


def find_circ_id(row):
    for candidate in ("circRNA", "circRNA_ID", "CircRNA", "ID", "row_id"):
        value = str(row.get(candidate, "")).strip()
        if value and looks_like_circ(value):
            return split_circ_gene(value)
    for value in row.values():
        if looks_like_circ(value):
            return split_circ_gene(value)
    return "", ""


def parse_bsj_id(circ_id):
    match = re.match(r"^([^:]+):(\d+)(?:\||-|\.\.)(\d+)$", str(circ_id).strip())
    if match is None:
        return None
    chrom, start_s, end_s = match.groups()
    start, end = int(start_s), int(end_s)
    if end < start:
        start, end = end, start
    return chrom, start, end


def read_ciri3_annotation(path):
    strand_map = {}
    gene_map = {}
    path = Path(path)
    if not path.exists():
        return strand_map, gene_map

    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            return strand_map, gene_map

        id_col = next(
            (col for col in ("circRNA_ID", "circRNA", "circ_id", "junction") if col in reader.fieldnames),
            None,
        )
        strand_col = next(
            (col for col in ("strand", "Strand", "bsj_strand") if col in reader.fieldnames),
            None,
        )
        gene_col = next(
            (col for col in ("gene_id", "gene", "gene_name") if col in reader.fieldnames),
            None,
        )
        if id_col is None:
            return strand_map, gene_map

        for row in reader:
            circ_id = str(row.get(id_col, "")).strip()
            if not circ_id:
                continue
            if strand_col is not None:
                strand = str(row.get(strand_col, "")).strip()
                if strand in {"+", "-"} and circ_id not in strand_map:
                    strand_map[circ_id] = strand
            if gene_col is not None:
                gene = str(row.get(gene_col, "")).strip()
                if gene and circ_id not in gene_map:
                    gene_map[circ_id] = gene
    return strand_map, gene_map


def result_metrics(row):
    _, pvalue_s = first_existing(
        row,
        ("PValue", "P.Value", "pvalue", "p_value", "Pvalue", "PVAL", "P.Val", "pval"),
    )
    _, padj_s = first_existing(
        row,
        ("padj", "FDR", "adj.P.Val", "qvalue", "q_value", "QValue", "BH"),
    )
    _, effect_s = first_existing(
        row,
        ("log2FoldChange", "logFC", "log2FC", "log2_fold_change", "Log2FC", "DeltaPSI"),
    )
    return parse_float(pvalue_s), parse_float(padj_s), parse_float(effect_s)


def row_to_event(row, method, strand_map, gene_map, padj_cutoff, lfc_cutoff):
    circ_id, gene_from_id = find_circ_id(row)
    parsed = parse_bsj_id(circ_id)
    if parsed is None:
        return None

    pvalue, padj, effect = result_metrics(row)
    if padj is None:
        return None
    if padj >= padj_cutoff:
        return None
    if method == "deseq2" and (effect is None or abs(effect) < lfc_cutoff):
        return None
    if method != "deseq2" and effect is not None and abs(effect) < lfc_cutoff:
        return None

    chrom, start, end = parsed
    _, row_strand = first_existing(row, ("strand", "Strand", "bsj_strand"))
    strand = strand_map.get(circ_id) or row_strand
    if strand not in {"+", "-"}:
        strand = "+"
    gene_id = gene_map.get(circ_id) or gene_from_id

    return {
        "circRNA_ID": circ_id,
        "gene_id": gene_id,
        "chrom": chrom,
        "start": start,
        "end": end,
        "strand": strand,
        "pvalue": pvalue,
        "padj": padj,
        "effect": effect,
    }


event_result = Path(snakemake.input.result)
ciri3_annotation = Path(snakemake.input.ciri3)
gff3 = Path(snakemake.input.gff3).resolve()
outdir = Path(snakemake.params.outdir)
plots_root = outdir / "plots"
manifest_path = Path(snakemake.output.manifest)
done_path = Path(snakemake.output.done)
log_path = Path(snakemake.log[0])

method = str(snakemake.wildcards.method)
comparison = str(snakemake.wildcards.comparison)
padj_cutoff = float(snakemake.params.padj_cutoff)
lfc_cutoff = float(snakemake.params.lfc_cutoff)
max_events = int(snakemake.params.max_events)
bsj_flank = int(snakemake.params.bsj_flank)
fail_on_error = bool(snakemake.params.fail_on_error)

if outdir.exists():
    shutil.rmtree(outdir)
plots_root.mkdir(parents=True, exist_ok=True)
manifest_path.parent.mkdir(parents=True, exist_ok=True)
done_path.parent.mkdir(parents=True, exist_ok=True)
log_path.parent.mkdir(parents=True, exist_ok=True)

_, rows = read_tabular_loose(event_result)
strand_map, gene_map = read_ciri3_annotation(ciri3_annotation)

events = []
for row in rows:
    event = row_to_event(row, method, strand_map, gene_map, padj_cutoff, lfc_cutoff)
    if event is not None:
        events.append(event)

events.sort(
    key=lambda event: (
        event["padj"] if event["padj"] is not None else float("inf"),
        event["pvalue"] if event["pvalue"] is not None else float("inf"),
        -(abs(event["effect"]) if event["effect"] is not None else 0),
        event["circRNA_ID"],
    )
)

deduped = []
seen = set()
for event in events:
    if event["circRNA_ID"] in seen:
        continue
    seen.add(event["circRNA_ID"])
    deduped.append(event)
events = deduped[:max_events] if max_events > 0 else deduped

base_cmd = [
    sys.executable,
    str(Path(snakemake.input.script).resolve()),
    "--b1",
    str(Path(snakemake.input.b1).resolve()),
    "--b2",
    str(Path(snakemake.input.b2).resolve()),
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
extra_args = str(snakemake.params.extra_args or "")
if extra_args:
    base_cmd.extend(shlex.split(extra_args))

env = os.environ.copy()
env.setdefault("MPLBACKEND", "Agg")
manifest_rows = []

with open(log_path, "w") as log_handle:
    log_handle.write(
        f"Input BSJ DEG result: {event_result}\n"
        f"Method: {method}\n"
        f"Comparison: {comparison}\n"
        f"Selected BSJs: {len(events)}\n"
        f"Thresholds: padj/FDR < {padj_cutoff}, |effect| >= {lfc_cutoff} "
        f"when an effect column is available; flank={bsj_flank}\n\n"
    )

    for index, event in enumerate(events, start=1):
        region_start = max(1, int(event["start"]) - bsj_flank)
        region_end = int(event["end"]) + bsj_flank
        coordinate = (
            f"{event['chrom']}:{event['strand']}:{region_start}:{region_end}:{gff3}"
        )
        event_key = (
            f"{index:05d}_{safe_name(event['gene_id'] or 'gene')}_"
            f"{safe_name(event['circRNA_ID'])}"
        )
        plot_outdir = plots_root / event_key
        plot_outdir.mkdir(parents=True, exist_ok=True)

        cmd = base_cmd + [
            "-c",
            coordinate,
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

        pdfs = sorted(str(path.resolve()) for path in plot_outdir.glob("**/*.pdf"))
        status = "ok" if proc.returncode == 0 and pdfs else "failed"
        if status == "failed" and fail_on_error:
            raise RuntimeError(
                f"rmats2sashimiplot failed for {event_key}; see {log_path}"
            )
        manifest_rows.append(
            {
                "method": method,
                "comparison": comparison,
                "rank": index,
                "circRNA_ID": event["circRNA_ID"],
                "gene_id": event["gene_id"],
                "chrom": event["chrom"],
                "bsj_start": event["start"],
                "bsj_end": event["end"],
                "strand": event["strand"],
                "region_start": region_start,
                "region_end": region_end,
                "pvalue": "" if event["pvalue"] is None else event["pvalue"],
                "padj": "" if event["padj"] is None else event["padj"],
                "effect": "" if event["effect"] is None else event["effect"],
                "coordinate": coordinate,
                "plot_dir": str(plot_outdir.resolve()),
                "pdfs": ",".join(pdfs),
                "status": status,
                "returncode": proc.returncode,
            }
        )
        log_handle.write(f"[{event_key}] status={status} pdfs={len(pdfs)} returncode={proc.returncode}\n\n")

fieldnames = [
    "method",
    "comparison",
    "rank",
    "circRNA_ID",
    "gene_id",
    "chrom",
    "bsj_start",
    "bsj_end",
    "strand",
    "region_start",
    "region_end",
    "pvalue",
    "padj",
    "effect",
    "coordinate",
    "plot_dir",
    "pdfs",
    "status",
    "returncode",
]
with open(manifest_path, "w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(manifest_rows)

done_path.write_text(
    f"events_selected\t{len(events)}\n"
    f"plots_ok\t{sum(1 for row in manifest_rows if row['status'] == 'ok')}\n"
    f"plots_failed\t{sum(1 for row in manifest_rows if row['status'] == 'failed')}\n"
)
