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

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pysam
from matplotlib.patches import Arc


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


def format_count(value):
    value = float(value)
    rounded = round(value)
    if abs(value - rounded) < 1e-6:
        return str(int(rounded))
    return f"{value:.2f}".rstrip("0").rstrip(".")


def format_ratio(value):
    return f"{float(value):.3f}".rstrip("0").rstrip(".")


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


def read_sashimi_samples(path):
    path = Path(path)
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Empty sashimi samples table: {path}")
        if "sample" not in reader.fieldnames or "group" not in reader.fieldnames:
            raise ValueError(f"Sashimi samples table must contain sample and group columns: {path}")

        groups = []
        group_to_samples = {}
        samples = []
        sample_to_bam = {}
        sample_to_index = {}
        for row in reader:
            sample = str(row.get("sample", "")).strip()
            group = str(row.get("group", "")).strip()
            if not sample or not group:
                continue
            if group not in group_to_samples:
                group_to_samples[group] = []
                groups.append(group)
            group_to_samples[group].append(sample)
            samples.append(sample)
            bam = str(row.get("bam", "")).strip()
            if bam:
                sample_to_bam[sample] = bam
            index = str(row.get("rmats2sashimi_index", "")).strip()
            if index:
                sample_to_index[sample] = index

    if not groups:
        raise ValueError(f"No sample/group records in sashimi samples table: {path}")
    return {
        "groups": groups,
        "group_to_samples": group_to_samples,
        "samples": samples,
        "sample_to_bam": sample_to_bam,
        "sample_to_index": sample_to_index,
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


def bsj_signal_summary(event):
    if "bsj_total_bsj" not in event:
        return ""
    return (
        f"CIRI3_BSJ={format_count(event['bsj_total_bsj'])}, "
        f"CIRI3_FSJ={format_count(event['bsj_total_fsj'])}, "
        f"JR={format_ratio(event['bsj_junction_ratio'])}, "
        f"matrix_row={event['bsj_matrix_row_found']}"
    )


def bsj_display_name(event, base_name):
    summary = bsj_signal_summary(event)
    if not summary:
        return base_name
    short_summary = (
        f"BSJ-{format_count(event['bsj_total_bsj'])}_"
        f"FSJ-{format_count(event['bsj_total_fsj'])}_"
        f"JR-{format_ratio(event['bsj_junction_ratio'])}"
    )
    return f"{base_name}_{short_summary}"


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
        signal_note = bsj_signal_summary(event)
        annotation_note = "Synthetic BSJ annotation from circRNA DEG result"
        if signal_note:
            annotation_note = f"{annotation_note}; {signal_note}"
        span_note = "Back-splice junction span"
        if signal_note:
            span_note = f"{span_note}; {signal_note}"
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
                    "Name": bsj_display_name(event, name),
                    "Alias": event["circRNA_ID"],
                    "Note": annotation_note,
                },
            ),
            (
                "mRNA",
                start,
                end,
                {
                    "ID": f"{feature_id}.tx",
                    "Parent": f"{feature_id}.gene",
                    "Name": bsj_display_name(event, event["circRNA_ID"]),
                    "Alias": name,
                    "Note": span_note,
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


def matrix_locus_keys(chrom, start, end):
    chrom = str(chrom)
    keys = [(chrom, int(start), int(end))]
    if chrom.startswith("chr") and len(chrom) > 3:
        keys.append((chrom[3:], int(start), int(end)))
    else:
        keys.append((f"chr{chrom}", int(start), int(end)))
    return keys


def read_ciri3_matrix(path):
    path = Path(path)
    rows = {}
    by_locus = {}
    samples = []

    with path.open() as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if header is None or len(header) < 2:
            return {"path": path, "samples": samples, "rows": rows, "by_locus": by_locus}
        samples = [sample.strip() for sample in header[1:]]

        for row in reader:
            if not row:
                continue
            circ_id = row[0].strip()
            if not circ_id:
                continue
            values = {}
            for sample, raw in zip(samples, row[1:]):
                parsed = parse_float(raw)
                values[sample] = 0.0 if parsed is None else parsed
            for sample in samples[len(row) - 1 :]:
                values.setdefault(sample, 0.0)

            aliases = [circ_id]
            if ";" in circ_id:
                aliases.append(circ_id.split(";", 1)[0])
            for alias in aliases:
                rows.setdefault(alias, values)
                parsed = parse_bsj_id(alias)
                if parsed is None:
                    continue
                chrom, start, end = parsed
                for key in matrix_locus_keys(chrom, start, end):
                    by_locus.setdefault(key, values)

    return {"path": path, "samples": samples, "rows": rows, "by_locus": by_locus}


def lookup_ciri3_matrix_values(matrix, event):
    circ_id = str(event["circRNA_ID"])
    candidates = [
        circ_id,
        f"{event['chrom']}:{event['start']}|{event['end']}",
        f"{event['chrom']}:{event['start']}-{event['end']}",
        f"{event['chrom']}:{event['start']}..{event['end']}",
    ]
    for candidate in candidates:
        if candidate in matrix["rows"]:
            return matrix["rows"][candidate]

    for key in matrix_locus_keys(event["chrom"], event["start"], event["end"]):
        if key in matrix["by_locus"]:
            return matrix["by_locus"][key]
    return {}


def junction_ratio(bsj, fsj):
    denom = 2.0 * bsj + fsj
    if denom <= 0:
        return 0.0
    return (2.0 * bsj) / denom


def collect_ciri3_bsj_signal(event, bsj_matrix, fsj_matrix, sample_info):
    bsj_values = lookup_ciri3_matrix_values(bsj_matrix, event)
    fsj_values = lookup_ciri3_matrix_values(fsj_matrix, event)
    group_rows = []
    sample_rows = []
    total_bsj = 0.0
    total_fsj = 0.0

    for group in sample_info["groups"]:
        group_samples = sample_info["group_to_samples"][group]
        group_bsj = 0.0
        group_fsj = 0.0
        positive_bsj = 0
        for sample in group_samples:
            bsj = float(bsj_values.get(sample, 0.0))
            fsj = float(fsj_values.get(sample, 0.0))
            group_bsj += bsj
            group_fsj += fsj
            if bsj > 0:
                positive_bsj += 1
            sample_rows.append(
                {
                    "group": group,
                    "sample": sample,
                    "bsj_count": bsj,
                    "fsj_count": fsj,
                    "junction_ratio": junction_ratio(bsj, fsj),
                }
            )

        total_bsj += group_bsj
        total_fsj += group_fsj
        group_rows.append(
            {
                "group": group,
                "samples": list(group_samples),
                "n_samples": len(group_samples),
                "n_bsj_positive": positive_bsj,
                "bsj_count": group_bsj,
                "fsj_count": group_fsj,
                "junction_ratio": junction_ratio(group_bsj, group_fsj),
            }
        )

    return {
        "group_rows": group_rows,
        "sample_rows": sample_rows,
        "total_bsj": total_bsj,
        "total_fsj": total_fsj,
        "junction_ratio": junction_ratio(total_bsj, total_fsj),
        "matrix_row_found": bool(bsj_values or fsj_values),
    }


def sample_signal_by_name(signal):
    return {row["sample"]: row for row in signal["sample_rows"]}


def event_anchor_length(event):
    start = int(event["start"])
    end = int(event["end"])
    span = max(1, end - start + 1)
    max_left = max(1, start - 1)
    return max(1, min(50, max_left, max(8, span // 10)))


def event_bam_header(source_bam, event, region_end, anchor):
    chrom = str(event["chrom"])
    minimum_length = max(int(region_end), int(event["end"]) + anchor + 100)
    header = None
    if source_bam:
        source_path = Path(source_bam)
        if source_path.exists():
            try:
                with pysam.AlignmentFile(str(source_path), "rb") as src:
                    header = src.header.to_dict()
            except (OSError, ValueError):
                header = None

    if header is None:
        header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": []}
    header.setdefault("HD", {}).setdefault("VN", "1.6")
    header.setdefault("SQ", [])

    sq_by_name = {entry.get("SN"): entry for entry in header["SQ"]}
    if chrom not in sq_by_name:
        header["SQ"].append({"SN": chrom, "LN": minimum_length})
    else:
        try:
            current_length = int(sq_by_name[chrom].get("LN", 0))
        except (TypeError, ValueError):
            current_length = 0
        if current_length < minimum_length:
            sq_by_name[chrom]["LN"] = minimum_length
    return header


def write_ciri3_pseudo_junction_reads(out_bam, event, bsj_count, read_prefix):
    count = max(0, int(round(float(bsj_count))))
    if count <= 0:
        return 0

    chrom = str(event["chrom"])
    start = int(event["start"])
    end = int(event["end"])
    if end - start <= 1:
        return 0

    anchor = event_anchor_length(event)
    skip = end - start - 1
    sequence = "A" * (anchor * 2)
    qualities = pysam.qualitystring_to_array("I" * len(sequence))
    reference_id = out_bam.get_tid(chrom)
    if reference_id < 0:
        return 0

    for read_index in range(1, count + 1):
        aln = pysam.AlignedSegment()
        aln.query_name = f"{read_prefix}_{read_index}"
        aln.query_sequence = sequence
        aln.flag = 16 if event.get("strand") == "-" else 0
        aln.reference_id = reference_id
        aln.reference_start = start - anchor
        aln.mapping_quality = 255
        aln.cigartuples = [(0, anchor), (3, skip), (0, anchor)]
        aln.query_qualities = qualities
        aln.set_tag("NH", 1, value_type="i")
        aln.set_tag("XS", event.get("strand", "+"), value_type="A")
        aln.set_tag("ZB", "CIRI3_BSJ", value_type="Z")
        out_bam.write(aln)
    return count


def write_event_bam(
    event,
    signal,
    sample,
    source_bam,
    out_bam,
    region_start,
    region_end,
    include_real_reads,
    log_handle,
):
    out_bam = Path(out_bam)
    out_bam.parent.mkdir(parents=True, exist_ok=True)
    unsorted_bam = out_bam.with_suffix(".unsorted.bam")
    anchor = event_anchor_length(event)
    header = event_bam_header(source_bam, event, region_end, anchor)
    chrom = str(event["chrom"])
    copied_reads = 0
    sample_signal = sample_signal_by_name(signal).get(sample, {})
    bsj_count = float(sample_signal.get("bsj_count", 0.0))

    with pysam.AlignmentFile(str(unsorted_bam), "wb", header=header) as out:
        if include_real_reads and source_bam:
            source_path = Path(source_bam)
            if source_path.exists():
                try:
                    with pysam.AlignmentFile(str(source_path), "rb") as src:
                        if chrom in src.references:
                            for read in src.fetch(chrom, max(0, int(region_start) - 1), int(region_end)):
                                out.write(read)
                                copied_reads += 1
                        else:
                            log_handle.write(
                                f"[{event['circRNA_ID']}][{sample}] source BAM lacks chrom {chrom}: "
                                f"{source_path}\n"
                            )
                except (OSError, ValueError) as exc:
                    log_handle.write(
                        f"[{event['circRNA_ID']}][{sample}] could not copy real BAM reads from "
                        f"{source_path}: {exc}\n"
                    )

        synthetic_reads = write_ciri3_pseudo_junction_reads(
            out,
            event,
            bsj_count,
            f"{safe_name(sample)}_{safe_name(event['circRNA_ID'])}_CIRI3BSJ",
        )

    pysam.sort("-o", str(out_bam), str(unsorted_bam))
    pysam.index(str(out_bam))
    unsorted_bam.unlink(missing_ok=True)
    return {
        "bam": str(out_bam.resolve()),
        "source_bam": str(source_bam or ""),
        "copied_reads": copied_reads,
        "synthetic_bsj_reads": synthetic_reads,
    }


def write_rmats2sashimi_inputs(input_dir, sample_info, sample_to_bam):
    input_dir = Path(input_dir)
    input_dir.mkdir(parents=True, exist_ok=True)
    group_names = list(sample_info["groups"])
    if len(group_names) < 2:
        raise ValueError("CIRI3 BSJ MISO input generation requires at least two groups.")

    first_group_samples = list(sample_info["group_to_samples"][group_names[0]])
    other_samples = [
        sample
        for group in group_names[1:]
        for sample in sample_info["group_to_samples"][group]
    ]
    ordered_samples = first_group_samples + other_samples
    sample_to_index = {sample: index for index, sample in enumerate(ordered_samples, start=1)}

    b1_path = input_dir / "b1.txt"
    b2_path = input_dir / "b2.txt"
    group_info_path = input_dir / "grouping.gf"
    samples_path = input_dir / "samples.tsv"

    b1_path.write_text(",".join(sample_to_bam[sample] for sample in first_group_samples) + "\n")
    b2_path.write_text(",".join(sample_to_bam[sample] for sample in other_samples) + "\n")
    with group_info_path.open("w") as handle:
        for group in group_names:
            indices = [str(sample_to_index[sample]) for sample in sample_info["group_to_samples"][group]]
            handle.write(f"{group}: {','.join(indices)}\n")
    with samples_path.open("w") as handle:
        handle.write("sample\tgroup\trmats2sashimi_index\tbam\n")
        for group in group_names:
            for sample in sample_info["group_to_samples"][group]:
                handle.write(
                    f"{sample}\t{group}\t{sample_to_index[sample]}\t{sample_to_bam[sample]}\n"
                )

    return {
        "b1": str(b1_path.resolve()),
        "b2": str(b2_path.resolve()),
        "group_info": str(group_info_path.resolve()),
        "samples": str(samples_path.resolve()),
    }


def build_ciri3_miso_inputs(
    event,
    signal,
    sample_info,
    event_key,
    plot_outdir,
    region_start,
    region_end,
    include_real_reads,
    log_handle,
):
    input_root = Path(plot_outdir) / "ciri3_miso_inputs"
    bam_root = input_root / "bams"
    sample_to_bam = {}
    bam_rows = []
    for index, sample in enumerate(sample_info["samples"], start=1):
        out_bam = bam_root / f"{index:03d}_{safe_name(sample)}.bam"
        source_bam = sample_info["sample_to_bam"].get(sample, "")
        row = write_event_bam(
            event,
            signal,
            sample,
            source_bam,
            out_bam,
            region_start,
            region_end,
            include_real_reads,
            log_handle,
        )
        sample_to_bam[sample] = row["bam"]
        bam_rows.append(row)

    input_files = write_rmats2sashimi_inputs(input_root, sample_info, sample_to_bam)
    synthetic_reads = sum(row["synthetic_bsj_reads"] for row in bam_rows)
    copied_reads = sum(row["copied_reads"] for row in bam_rows)
    log_handle.write(
        f"[{event_key}][ciri3_miso_inputs] mode="
        f"{'real_plus_ciri3_bsj' if include_real_reads else 'ciri3_bsj_only'} "
        f"synthetic_bsj_reads={synthetic_reads} copied_real_reads={copied_reads} "
        f"input_dir={input_root.resolve()}\n"
    )
    return {
        **input_files,
        "bam_source": "real_BAM_plus_CIRI3_BSJ_pseudo_reads"
        if include_real_reads
        else "CIRI3_BSJ_pseudo_reads",
        "synthetic_bsj_reads": synthetic_reads,
        "copied_real_reads": copied_reads,
        "input_dir": str(input_root.resolve()),
    }


def write_ciri3_signal_table(path, signal):
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "level",
        "group",
        "sample",
        "n_samples",
        "n_bsj_positive",
        "bsj_count",
        "fsj_count",
        "junction_ratio",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in signal["group_rows"]:
            writer.writerow(
                {
                    "level": "group",
                    "group": row["group"],
                    "sample": "",
                    "n_samples": row["n_samples"],
                    "n_bsj_positive": row["n_bsj_positive"],
                    "bsj_count": format_count(row["bsj_count"]),
                    "fsj_count": format_count(row["fsj_count"]),
                    "junction_ratio": format_ratio(row["junction_ratio"]),
                }
            )
        for row in signal["sample_rows"]:
            writer.writerow(
                {
                    "level": "sample",
                    "group": row["group"],
                    "sample": row["sample"],
                    "n_samples": "",
                    "n_bsj_positive": "",
                    "bsj_count": format_count(row["bsj_count"]),
                    "fsj_count": format_count(row["fsj_count"]),
                    "junction_ratio": format_ratio(row["junction_ratio"]),
                }
            )


def run_ciri3_bsj_signal_plot(
    event,
    signal,
    region_start,
    region_end,
    plot_outdir,
    fig_width,
    fig_height,
    font_size,
    colors,
    event_key,
    log_handle,
    plot_label,
):
    try:
        plot_outdir.mkdir(parents=True, exist_ok=True)
        sashimi_dir = plot_outdir / "Sashimi_plot"
        sashimi_dir.mkdir(parents=True, exist_ok=True)
        source_data = plot_outdir / "ciri3_bsj_fsj_signal.tsv"
        write_ciri3_signal_table(source_data, signal)

        group_rows = signal["group_rows"]
        group_count = max(1, len(group_rows))
        width = max(float(fig_width), 8.0)
        height = max(float(fig_height) if fig_height > 0 else 0.0, 2.8 + group_count * 0.85)
        height = min(max(height, 4.0), 14.0)
        font_size = max(int(font_size), 6)

        fig, ax = plt.subplots(figsize=(width, height))
        start = int(event["start"])
        end = int(event["end"])
        if end <= start:
            end = start + 1
        x_min = min(int(region_start), start)
        x_max = max(int(region_end), end)
        if x_max <= x_min:
            x_max = x_min + 1
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(0, group_count + 1.2)

        max_bsj = max((row["bsj_count"] for row in group_rows), default=0.0)
        max_fsj = max((row["fsj_count"] for row in group_rows), default=0.0)
        label_font = max(font_size - 1, 6)

        for idx, row in enumerate(group_rows):
            y = group_count - idx
            color = colors[idx % len(colors)] if colors else f"C{idx}"
            fsj = row["fsj_count"]
            bsj = row["bsj_count"]
            if fsj > 0:
                fsj_width = 0.8 + 3.0 * math.sqrt(fsj / max_fsj) if max_fsj > 0 else 0.8
                ax.plot(
                    [start, end],
                    [y, y],
                    color=color,
                    linewidth=fsj_width,
                    alpha=0.28,
                    solid_capstyle="round",
                )
            else:
                ax.plot([start, end], [y, y], color="0.82", linewidth=0.8, alpha=0.7)

            if bsj > 0:
                arc_width = 1.2 + 5.0 * math.sqrt(bsj / max_bsj) if max_bsj > 0 else 1.2
                arc_height = max(0.45, 0.18 + min(0.75, math.log2(bsj + 1) * 0.18))
                arc = Arc(
                    ((start + end) / 2.0, y),
                    width=(end - start),
                    height=arc_height,
                    theta1=0,
                    theta2=180,
                    linewidth=arc_width,
                    color=color,
                    alpha=0.95,
                )
                ax.add_patch(arc)
                ax.text(
                    (start + end) / 2.0,
                    y + arc_height / 2.0 + 0.08,
                    format_count(bsj),
                    ha="center",
                    va="bottom",
                    fontsize=label_font,
                    color=color,
                )

            ax.vlines([start, end], y - 0.12, y + 0.12, color="0.25", linewidth=0.8)
            ax.text(
                -0.02,
                y,
                row["group"],
                transform=ax.get_yaxis_transform(),
                ha="right",
                va="center",
                fontsize=font_size,
            )
            ax.text(
                1.01,
                y,
                (
                    f"BSJ={format_count(bsj)}  FSJ={format_count(fsj)}  "
                    f"JR={format_ratio(row['junction_ratio'])}  "
                    f"n+={row['n_bsj_positive']}/{row['n_samples']}"
                ),
                transform=ax.get_yaxis_transform(),
                ha="left",
                va="center",
                fontsize=label_font,
            )

        ax.axvline(start, color="0.35", linewidth=0.8, linestyle=":")
        ax.axvline(end, color="0.35", linewidth=0.8, linestyle=":")
        ax.set_yticks([])
        ax.set_xlabel(f"{event['chrom']} coordinate", fontsize=label_font)
        ax.set_title(
            (
                f"{event['gene_id'] or event['circRNA_ID']} | {event['circRNA_ID']}\n"
                "CIRI3 BSJ_Matrix / FSJ_Matrix signal"
            ),
            fontsize=font_size + 1,
        )
        ax.text(
            0.01,
            0.02,
            (
                "Arc: CIRI3 BSJ count; baseline: CIRI3 FSJ count; "
                "JR=2*BSJ/(2*BSJ+FSJ)"
            ),
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=label_font,
            color="0.25",
        )
        if not signal["matrix_row_found"]:
            ax.text(
                0.5,
                0.5,
                "No matching CIRI3 BSJ/FSJ matrix row",
                transform=ax.transAxes,
                ha="center",
                va="center",
                fontsize=font_size,
                color="0.35",
            )
        for spine in ["left", "right", "top"]:
            ax.spines[spine].set_visible(False)
        fig.tight_layout()

        pdf_path = sashimi_dir / f"{event_key}.ciri3_bsj_fsj_signal.pdf"
        fig.savefig(pdf_path, bbox_inches="tight")
        plt.close(fig)
        pdfs = [str(pdf_path.resolve())] if pdf_path.is_file() else []
        status = "ok" if pdfs else "failed"
        log_handle.write(
            f"[{event_key}][{plot_label}_bsj_signal] status={status} "
            f"pdfs={len(pdfs)} total_bsj={format_count(signal['total_bsj'])} "
            f"total_fsj={format_count(signal['total_fsj'])}\n\n"
        )
        return {
            "plot_dir": str(plot_outdir.resolve()),
            "pdfs": pdfs,
            "source_data": str(source_data.resolve()),
            "status": status,
            "returncode": 0 if status == "ok" else 1,
        }
    except Exception as exc:
        log_handle.write(f"[{event_key}][{plot_label}_bsj_signal] failed: {exc}\n\n")
        return {
            "plot_dir": str(plot_outdir.resolve()),
            "pdfs": [],
            "source_data": "",
            "status": "failed",
            "returncode": 1,
        }


def combine_plot_results(miso_result, signal_result):
    combined = dict(miso_result)
    combined["pdfs"] = sorted(set(miso_result["pdfs"] + signal_result["pdfs"]))
    combined["signal_pdfs"] = signal_result["pdfs"]
    combined["signal_source_data"] = signal_result["source_data"]
    combined["status"] = "ok" if miso_result["status"] == "ok" and signal_result["status"] == "ok" else "failed"
    combined["returncode"] = max(int(miso_result["returncode"]), int(signal_result["returncode"]))
    return combined


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
bsj_matrix_path = Path(snakemake.input.bsj_matrix)
fsj_matrix_path = Path(snakemake.input.fsj_matrix)
samples_path = Path(snakemake.input.samples)
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
sample_info = read_sashimi_samples(samples_path)
ciri3_bsj_matrix = read_ciri3_matrix(bsj_matrix_path)
ciri3_fsj_matrix = read_ciri3_matrix(fsj_matrix_path)

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
for event in events:
    signal = collect_ciri3_bsj_signal(
        event,
        ciri3_bsj_matrix,
        ciri3_fsj_matrix,
        sample_info,
    )
    event["bsj_signal"] = signal
    event["bsj_total_bsj"] = signal["total_bsj"]
    event["bsj_total_fsj"] = signal["total_fsj"]
    event["bsj_junction_ratio"] = signal["junction_ratio"]
    event["bsj_matrix_row_found"] = "yes" if signal["matrix_row_found"] else "no"
write_augmented_gff3(base_gff3, augmented_gff3, events)

colors = list(snakemake.params.colors)
plot_colors = []
if colors:
    plot_colors = colors[: group_summary["groups"]]
extra_args = str(snakemake.params.extra_args or "")


def make_base_cmd(b1_path, b2_path, group_info_path):
    cmd = [
        sys.executable,
        str(Path(snakemake.input.script).resolve()),
        "--b1",
        str(Path(b1_path).resolve()),
        "--b2",
        str(Path(b2_path).resolve()),
        "--l1",
        str(snakemake.params.label_b1),
        "--l2",
        str(snakemake.params.label_b2),
        "--exon_s",
        str(snakemake.params.exon_s),
        "--intron_s",
        str(snakemake.params.intron_s),
        "--group-info",
        str(Path(group_info_path).resolve()),
        "--min-counts",
        str(snakemake.params.min_counts),
    ]
    if colors and len(plot_colors) == group_summary["groups"]:
        cmd.extend(["--color", ",".join(plot_colors)])
    if extra_args:
        cmd.extend(shlex.split(extra_args))
    return cmd

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
        f"BSJ labels: CIRI3 BSJ_Matrix / FSJ_Matrix fields embedded in synthetic GFF3\n"
        f"CIRI3 BSJ loading: event-level pseudo-BAM junction reads are generated "
        f"and passed to rmats2sashimiplot/MISO\n"
        f"CIRI3 BSJ matrix: {bsj_matrix_path.resolve()}\n"
        f"CIRI3 FSJ matrix: {fsj_matrix_path.resolve()}\n"
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
        bsj_signal = event["bsj_signal"]
        full_inputs = build_ciri3_miso_inputs(
            event,
            bsj_signal,
            sample_info,
            event_key,
            plots_root / event_key,
            region_start,
            region_end,
            True,
            log_handle,
        )
        full_miso_plot = run_sashimi_plot(
            make_base_cmd(full_inputs["b1"], full_inputs["b2"], full_inputs["group_info"]),
            style_args,
            coordinate,
            plots_root / event_key,
            log_handle,
            env,
            event_key,
            "full",
        )
        full_signal_plot = run_ciri3_bsj_signal_plot(
            event,
            bsj_signal,
            region_start,
            region_end,
            plots_root / event_key,
            fig_width,
            fig_height,
            font_size,
            colors,
            event_key,
            log_handle,
            "full",
        )
        full_plot = combine_plot_results(full_miso_plot, full_signal_plot)
        bsj_only_inputs = build_ciri3_miso_inputs(
            event,
            bsj_signal,
            sample_info,
            event_key,
            bsj_only_plots_root / event_key,
            region_start,
            region_end,
            False,
            log_handle,
        )
        bsj_only_miso_plot = run_sashimi_plot(
            make_base_cmd(
                bsj_only_inputs["b1"],
                bsj_only_inputs["b2"],
                bsj_only_inputs["group_info"],
            ),
            style_args,
            bsj_only_coordinate,
            bsj_only_plots_root / event_key,
            log_handle,
            env,
            event_key,
            "bsj_only",
        )
        bsj_only_signal_plot = run_ciri3_bsj_signal_plot(
            event,
            bsj_signal,
            region_start,
            region_end,
            bsj_only_plots_root / event_key,
            fig_width,
            fig_height,
            font_size,
            colors,
            event_key,
            log_handle,
            "bsj_only",
        )
        bsj_only_plot = combine_plot_results(bsj_only_miso_plot, bsj_only_signal_plot)
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
                "bsj_total_bsj": format_count(bsj_signal["total_bsj"]),
                "bsj_total_fsj": format_count(bsj_signal["total_fsj"]),
                "bsj_junction_ratio": format_ratio(bsj_signal["junction_ratio"]),
                "bsj_matrix_row_found": event["bsj_matrix_row_found"],
                "bsj_signal_pdfs": ",".join(full_plot["signal_pdfs"]),
                "bsj_signal_source_data": full_plot["signal_source_data"],
                "miso_bam_source": full_inputs["bam_source"],
                "miso_input_dir": full_inputs["input_dir"],
                "miso_input_b1": full_inputs["b1"],
                "miso_input_b2": full_inputs["b2"],
                "miso_input_group_info": full_inputs["group_info"],
                "miso_synthetic_bsj_reads": full_inputs["synthetic_bsj_reads"],
                "miso_copied_real_reads": full_inputs["copied_real_reads"],
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
                "bsj_only_source": "rmats2sashimiplot_MISO_from_CIRI3_BSJ_pseudo_BAM",
                "bsj_only_bsj_matrix": str(bsj_matrix_path.resolve()),
                "bsj_only_fsj_matrix": str(fsj_matrix_path.resolve()),
                "bsj_only_miso_bam_source": bsj_only_inputs["bam_source"],
                "bsj_only_miso_input_dir": bsj_only_inputs["input_dir"],
                "bsj_only_miso_input_b1": bsj_only_inputs["b1"],
                "bsj_only_miso_input_b2": bsj_only_inputs["b2"],
                "bsj_only_miso_input_group_info": bsj_only_inputs["group_info"],
                "bsj_only_miso_synthetic_bsj_reads": bsj_only_inputs["synthetic_bsj_reads"],
                "bsj_only_miso_copied_real_reads": bsj_only_inputs["copied_real_reads"],
                "bsj_only_source_data": bsj_only_plot["signal_source_data"],
                "bsj_only_signal_pdfs": ",".join(bsj_only_plot["signal_pdfs"]),
                "bsj_only_total_bsj": format_count(bsj_signal["total_bsj"]),
                "bsj_only_total_fsj": format_count(bsj_signal["total_fsj"]),
                "bsj_only_junction_ratio": format_ratio(bsj_signal["junction_ratio"]),
                "bsj_only_matrix_row_found": event["bsj_matrix_row_found"],
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
    "bsj_total_bsj",
    "bsj_total_fsj",
    "bsj_junction_ratio",
    "bsj_matrix_row_found",
    "bsj_signal_pdfs",
    "bsj_signal_source_data",
    "miso_bam_source",
    "miso_input_dir",
    "miso_input_b1",
    "miso_input_b2",
    "miso_input_group_info",
    "miso_synthetic_bsj_reads",
    "miso_copied_real_reads",
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
    "bsj_only_source",
    "bsj_only_bsj_matrix",
    "bsj_only_fsj_matrix",
    "bsj_only_miso_bam_source",
    "bsj_only_miso_input_dir",
    "bsj_only_miso_input_b1",
    "bsj_only_miso_input_b2",
    "bsj_only_miso_input_group_info",
    "bsj_only_miso_synthetic_bsj_reads",
    "bsj_only_miso_copied_real_reads",
    "bsj_only_source_data",
    "bsj_only_signal_pdfs",
    "bsj_only_total_bsj",
    "bsj_only_total_fsj",
    "bsj_only_junction_ratio",
    "bsj_only_matrix_row_found",
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
