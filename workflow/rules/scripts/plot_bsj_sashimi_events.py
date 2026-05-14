import csv
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
from urllib.parse import quote


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


def clamp(value, minimum, maximum):
    return max(minimum, min(maximum, value))


def format_float(value):
    return f"{float(value):.2f}".rstrip("0").rstrip(".")


def snakemake_named_path(collection, name, fallback):
    value = getattr(collection, name, None)
    return Path(value) if value not in (None, "") else fallback


def gff3_attrs(attrs):
    return ";".join(
        f"{quote(str(key), safe='._:-')}={quote(str(value), safe='._:-')}"
        for key, value in attrs.items()
        if value not in (None, "")
    )


def parse_group_info(path):
    group_count = 0
    sample_indices = set()
    max_label_len = 0
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line or ":" not in line:
            continue
        group_name, indices_s = line.split(":", 1)
        group_count += 1
        max_label_len = max(max_label_len, len(group_name.strip()))
        for value in indices_s.split(","):
            value = value.strip()
            if value:
                sample_indices.add(value)
    return {
        "groups": max(group_count, 2),
        "samples": max(len(sample_indices), 2),
        "max_label_len": max_label_len,
    }


def scaled_plot_style(
    event,
    region_start,
    region_end,
    group_summary,
    auto_scale,
    base_width,
    base_height,
    min_width,
    max_width,
    min_height,
    max_height,
    base_font_size,
    min_font_size,
):
    if not auto_scale:
        return base_width, base_height, base_font_size

    span = max(1, int(region_end) - int(region_start) + 1)
    group_count = group_summary["groups"]
    sample_count = group_summary["samples"]
    label_len = max(
        len(str(event.get("circRNA_ID", ""))),
        len(str(event.get("gene_id", ""))),
        group_summary["max_label_len"],
    )

    width_by_span = 6.5 + math.log2(max(2.0, span / 100.0)) * 1.25
    width_by_label = min_width + max(0, label_len - 24) * 0.05
    width_by_groups = min_width + max(0, group_count - 2) * 0.35
    fig_width = clamp(
        max(base_width, min_width, width_by_span, width_by_label, width_by_groups),
        min_width,
        max_width,
    )

    height_by_groups = 3.0 + group_count * 1.15 + max(0, sample_count - group_count) * 0.12
    fig_height = clamp(
        max(base_height if base_height > 0 else 0, min_height, height_by_groups),
        min_height,
        max_height,
    )

    font_size = int(base_font_size)
    sample_density = sample_count / max(fig_width, 1)
    if group_count >= 4 or sample_count >= 8 or span < 1000 or sample_density > 0.8:
        font_size -= 1
    if group_count >= 6 or sample_count >= 14 or span < 400 or label_len >= 45:
        font_size -= 1
    font_size = max(int(min_font_size), font_size)

    return round(fig_width, 2), round(fig_height, 2), font_size


def synthetic_bsj_gff3_lines(events):
    lines = []
    for index, event in enumerate(events, start=1):
        start = int(event["start"])
        end = int(event["end"])
        span = max(1, end - start + 1)
        anchor = max(1, min(50, span // 10 if span >= 10 else span))
        left_end = min(end, start + anchor - 1)
        right_start = max(start, end - anchor + 1)
        feature_id = f"bsj-{index:05d}-{safe_name(event['circRNA_ID'])}"
        name = event["gene_id"] or event["circRNA_ID"]
        common = [
            str(event["chrom"]),
            "circRNA_smk",
            None,
            None,
            None,
            ".",
            str(event["strand"]),
            ".",
            None,
        ]
        rows = [
            (
                "gene",
                start,
                end,
                {
                    "ID": f"{feature_id}.gene",
                    "Name": name,
                    "Alias": event["circRNA_ID"],
                    "Note": "Synthetic BSJ annotation from circRNA DEG result",
                },
            ),
            (
                "mRNA",
                start,
                end,
                {
                    "ID": f"{feature_id}.tx",
                    "Parent": f"{feature_id}.gene",
                    "Name": event["circRNA_ID"],
                    "Alias": name,
                    "Note": "Back-splice junction span",
                },
            ),
            (
                "exon",
                start,
                left_end,
                {
                    "ID": f"{feature_id}.left_anchor",
                    "Parent": f"{feature_id}.tx",
                    "Name": "BSJ_left_anchor",
                },
            ),
            (
                "exon",
                right_start,
                end,
                {
                    "ID": f"{feature_id}.right_anchor",
                    "Parent": f"{feature_id}.tx",
                    "Name": "BSJ_right_anchor",
                },
            ),
        ]
        for feature, row_start, row_end, attrs in rows:
            parts = list(common)
            parts[2] = feature
            parts[3] = str(row_start)
            parts[4] = str(row_end)
            parts[8] = gff3_attrs(attrs)
            lines.append("\t".join(parts) + "\n")
    return lines


def write_augmented_gff3(base_gff3, augmented_gff3, events):
    synthetic_lines = synthetic_bsj_gff3_lines(events)
    augmented_gff3.parent.mkdir(parents=True, exist_ok=True)
    wrote_synthetic = False
    last_was_newline = True

    with Path(base_gff3).open() as src, augmented_gff3.open("w") as out:
        for line in src:
            if line.strip() == "##FASTA" and not wrote_synthetic:
                if not last_was_newline:
                    out.write("\n")
                out.write("###\n")
                out.write("# Synthetic BSJ annotations added by circRNA_smk.\n")
                out.writelines(synthetic_lines)
                out.write("###\n")
                wrote_synthetic = True
            out.write(line)
            last_was_newline = line.endswith("\n")
        if not wrote_synthetic:
            if not last_was_newline:
                out.write("\n")
            out.write("###\n")
            out.write("# Synthetic BSJ annotations added by circRNA_smk.\n")
            out.writelines(synthetic_lines)
            out.write("###\n")


def write_bsj_only_gff3(bsj_only_gff3, events):
    synthetic_lines = synthetic_bsj_gff3_lines(events)
    bsj_only_gff3.parent.mkdir(parents=True, exist_ok=True)
    with bsj_only_gff3.open("w") as out:
        out.write("##gff-version 3\n")
        out.write("# Synthetic BSJ-only annotations added by circRNA_smk.\n")
        if synthetic_lines:
            out.writelines(synthetic_lines)
            out.write("###\n")


def run_sashimi_plot(base_cmd, style_args, coordinate, plot_outdir, log_handle, env, event_key, plot_label):
    plot_outdir.mkdir(parents=True, exist_ok=True)
    cmd = base_cmd + style_args + [
        "-c",
        coordinate,
        "-o",
        str(plot_outdir.resolve()),
    ]
    log_handle.write(
        f"[{event_key}][{plot_label}] command: "
        f"{' '.join(shlex.quote(part) for part in cmd)}\n"
    )
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
    log_handle.write(
        f"[{event_key}][{plot_label}] status={status} "
        f"pdfs={len(pdfs)} returncode={proc.returncode}\n\n"
    )
    return {
        "coordinate": coordinate,
        "plot_dir": str(plot_outdir.resolve()),
        "pdfs": pdfs,
        "status": status,
        "returncode": proc.returncode,
    }


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
base_gff3 = Path(snakemake.input.gff3).resolve()
outdir = Path(snakemake.params.outdir)
plots_root = snakemake_named_path(snakemake.output, "plots", outdir / "plots")
bsj_only_plots_root = snakemake_named_path(
    snakemake.output,
    "bsj_only_plots",
    outdir / "plots_bsj_only",
)
annotation_root = outdir / "annotations"
bsj_only_annotation_root = annotation_root / "bsj_only"
augmented_gff3 = (annotation_root / "bsj_augmented.gff3").resolve()
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
auto_scale = bool(snakemake.params.auto_scale)
base_fig_width = float(snakemake.params.fig_width)
base_fig_height = float(snakemake.params.fig_height)
base_font_size = int(snakemake.params.font_size)
min_fig_width = float(snakemake.params.min_fig_width)
max_fig_width = float(snakemake.params.max_fig_width)
min_fig_height = float(snakemake.params.min_fig_height)
max_fig_height = float(snakemake.params.max_fig_height)
min_font_size = int(snakemake.params.min_font_size)

if outdir.exists():
    shutil.rmtree(outdir)
plots_root.mkdir(parents=True, exist_ok=True)
bsj_only_plots_root.mkdir(parents=True, exist_ok=True)
annotation_root.mkdir(parents=True, exist_ok=True)
bsj_only_annotation_root.mkdir(parents=True, exist_ok=True)
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
group_summary = parse_group_info(snakemake.input.group_info)
events = deduped[:max_events] if max_events > 0 else deduped
write_augmented_gff3(base_gff3, augmented_gff3, events)

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
]
colors = list(snakemake.params.colors)
plot_colors = []
if colors:
    plot_colors = colors[: group_summary["groups"]]
    if len(plot_colors) == group_summary["groups"]:
        base_cmd.extend(["--color", ",".join(plot_colors)])
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
        f"Base GFF3: {base_gff3}\n"
        f"Augmented GFF3: {augmented_gff3}\n"
        f"Full plot directory: {plots_root.resolve()}\n"
        f"BSJ-only plot directory: {bsj_only_plots_root.resolve()}\n"
        f"BSJ-only GFF3 directory: {bsj_only_annotation_root.resolve()}\n"
        f"Auto scale: {auto_scale}; groups={group_summary['groups']}; "
        f"samples={group_summary['samples']}\n"
        f"Colors passed to rmats2sashimiplot: {','.join(plot_colors) or 'default'}\n\n"
    )

    for index, event in enumerate(events, start=1):
        region_start = max(1, int(event["start"]) - bsj_flank)
        region_end = int(event["end"]) + bsj_flank
        fig_width, fig_height, font_size = scaled_plot_style(
            event,
            region_start,
            region_end,
            group_summary,
            auto_scale,
            base_fig_width,
            base_fig_height,
            min_fig_width,
            max_fig_width,
            min_fig_height,
            max_fig_height,
            base_font_size,
            min_font_size,
        )
        event_key = (
            f"{index:05d}_{safe_name(event['gene_id'] or 'gene')}_"
            f"{safe_name(event['circRNA_ID'])}"
        )
        coordinate = (
            f"{event['chrom']}:{event['strand']}:{region_start}:{region_end}:{augmented_gff3}"
        )
        bsj_only_gff3 = (bsj_only_annotation_root / f"{event_key}.gff3").resolve()
        write_bsj_only_gff3(bsj_only_gff3, [event])
        bsj_only_coordinate = (
            f"{event['chrom']}:{event['strand']}:{region_start}:{region_end}:{bsj_only_gff3}"
        )

        style_args = [
            "--font-size",
            str(font_size),
            "--fig-width",
            format_float(fig_width),
        ]
        if fig_height > 0:
            style_args.extend(["--fig-height", format_float(fig_height)])
        full_plot = run_sashimi_plot(
            base_cmd,
            style_args,
            coordinate,
            plots_root / event_key,
            log_handle,
            env,
            event_key,
            "full",
        )
        bsj_only_plot = run_sashimi_plot(
            base_cmd,
            style_args,
            bsj_only_coordinate,
            bsj_only_plots_root / event_key,
            log_handle,
            env,
            event_key,
            "bsj_only",
        )
        failed_variants = [
            label
            for label, result in (("full", full_plot), ("bsj_only", bsj_only_plot))
            if result["status"] == "failed"
        ]
        if failed_variants and fail_on_error:
            raise RuntimeError(
                f"rmats2sashimiplot failed for {event_key} "
                f"({', '.join(failed_variants)}); see {log_path}"
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
                "annotation_gff3": str(augmented_gff3),
                "bsj_only_coordinate": bsj_only_coordinate,
                "bsj_only_annotation_gff3": str(bsj_only_gff3),
                "fig_width": format_float(fig_width),
                "fig_height": format_float(fig_height) if fig_height > 0 else "",
                "font_size": font_size,
                "plot_dir": full_plot["plot_dir"],
                "pdfs": ",".join(full_plot["pdfs"]),
                "status": full_plot["status"],
                "returncode": full_plot["returncode"],
                "bsj_only_plot_dir": bsj_only_plot["plot_dir"],
                "bsj_only_pdfs": ",".join(bsj_only_plot["pdfs"]),
                "bsj_only_status": bsj_only_plot["status"],
                "bsj_only_returncode": bsj_only_plot["returncode"],
            }
        )

plots_root.mkdir(parents=True, exist_ok=True)
bsj_only_plots_root.mkdir(parents=True, exist_ok=True)

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
    "annotation_gff3",
    "bsj_only_coordinate",
    "bsj_only_annotation_gff3",
    "fig_width",
    "fig_height",
    "font_size",
    "plot_dir",
    "pdfs",
    "status",
    "returncode",
    "bsj_only_plot_dir",
    "bsj_only_pdfs",
    "bsj_only_status",
    "bsj_only_returncode",
]
with open(manifest_path, "w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(manifest_rows)

done_path.write_text(
    f"events_selected\t{len(events)}\n"
    f"plots_ok\t{sum(1 for row in manifest_rows if row['status'] == 'ok')}\n"
    f"plots_failed\t{sum(1 for row in manifest_rows if row['status'] == 'failed')}\n"
    f"bsj_only_plots_ok\t{sum(1 for row in manifest_rows if row['bsj_only_status'] == 'ok')}\n"
    f"bsj_only_plots_failed\t{sum(1 for row in manifest_rows if row['bsj_only_status'] == 'failed')}\n"
)
