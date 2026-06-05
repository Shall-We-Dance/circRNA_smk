import csv
import re
from collections import defaultdict
from pathlib import Path


samples = list(snakemake.params.samples)
circ_tables = list(map(Path, snakemake.input.circ_tables))
aggregate_bed = Path(snakemake.output.aggregate)
metadata_tsv = Path(snakemake.output.metadata)
per_sample_beds = list(map(Path, snakemake.output.per_sample))
log_path = Path(snakemake.log[0])


def safe_float(value, default=0.0):
    try:
        return float(str(value).strip())
    except (TypeError, ValueError):
        return default


def format_count(value):
    value = safe_float(value)
    if value.is_integer():
        return str(int(value))
    return f"{value:.6g}"


def bed_name(*parts):
    raw = "|".join(str(part) for part in parts if str(part) not in {"", "NA"})
    return re.sub(r"\s+", "_", raw) or "BSJ"


def bed_color(strand):
    if strand == "+":
        return "213,94,0"
    if strand == "-":
        return "0,114,178"
    return "80,80,80"


def bed_score(value, max_value):
    value = safe_float(value)
    if value <= 0 or max_value <= 0:
        return 0
    return max(1, min(1000, int(round(value / max_value * 1000))))


def bed12_row(row, count, max_count):
    start_1based = int(row["bsj_start"])
    end_1based = int(row["bsj_end"])
    if end_1based < start_1based:
        start_1based, end_1based = end_1based, start_1based

    chrom_start = max(0, start_1based - 1)
    chrom_end = max(chrom_start + 1, end_1based)
    score = bed_score(count, max_count)
    name = bed_name(
        row["circRNA_ID"],
        row.get("gene_id", "NA"),
        f"BSJ={format_count(count)}",
    )

    if chrom_end - chrom_start <= 1:
        block_count = 1
        block_sizes = "1"
        block_starts = "0"
    else:
        block_count = 2
        block_sizes = "1,1"
        block_starts = f"0,{chrom_end - chrom_start - 1}"

    return [
        row["chrom"],
        str(chrom_start),
        str(chrom_end),
        name,
        str(score),
        row.get("strand", ".") or ".",
        str(chrom_start),
        str(chrom_end),
        bed_color(row.get("strand", ".")),
        str(block_count),
        block_sizes,
        block_starts,
    ]


def read_circ_tables(paths):
    rows_by_sample = defaultdict(list)
    aggregate = {}
    sample_counts = defaultdict(dict)
    sample_fsj_counts = defaultdict(dict)

    for path in paths:
        with path.open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                circ_id = row.get("circRNA_ID", "").strip()
                if not circ_id:
                    continue
                sample = row.get("sample", "").strip()
                bsj_count = safe_float(row.get("bsj_count", 0))
                fsj_count = safe_float(row.get("fsj_count", 0))
                if sample:
                    rows_by_sample[sample].append(row)
                    sample_counts[circ_id][sample] = bsj_count
                    sample_fsj_counts[circ_id][sample] = fsj_count

                if circ_id not in aggregate:
                    aggregate[circ_id] = {
                        "row": row,
                        "total_bsj": 0.0,
                        "total_fsj": 0.0,
                        "samples_detected": 0,
                        "max_sample_bsj": 0.0,
                    }
                aggregate[circ_id]["total_bsj"] += bsj_count
                aggregate[circ_id]["total_fsj"] += fsj_count
                if bsj_count > 0:
                    aggregate[circ_id]["samples_detected"] += 1
                    aggregate[circ_id]["max_sample_bsj"] = max(
                        aggregate[circ_id]["max_sample_bsj"],
                        bsj_count,
                    )

    return rows_by_sample, aggregate, sample_counts, sample_fsj_counts


def write_bed(path, rows, label, count_key, max_count):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        handle.write(
            'track name="{}" description="CIRI3 BSJ junctions; '
            'BED score is scaled, raw counts are in all_samples.bsj.tsv" '
            'visibility=2 itemRgb="On" useScore=1\n'.format(label)
        )
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for row in rows:
            count = row[count_key]
            if safe_float(count) <= 0:
                continue
            writer.writerow(bed12_row(row, count, max_count))


def sample_count_string(counts, sample_names):
    parts = []
    for sample in sample_names:
        value = counts.get(sample, 0.0)
        parts.append(f"{sample}:{format_count(value)}")
    return ";".join(parts)


rows_by_sample, aggregate, sample_counts, sample_fsj_counts = read_circ_tables(circ_tables)

aggregate_rows = []
for circ_id, data in aggregate.items():
    row = dict(data["row"])
    row["total_bsj"] = data["total_bsj"]
    row["total_fsj"] = data["total_fsj"]
    row["samples_detected"] = data["samples_detected"]
    row["max_sample_bsj"] = data["max_sample_bsj"]
    aggregate_rows.append(row)

aggregate_rows.sort(
    key=lambda row: (row["chrom"], int(row["bsj_start"]), int(row["bsj_end"]), row["circRNA_ID"])
)
aggregate_max = max((safe_float(row["total_bsj"]) for row in aggregate_rows), default=0.0)
write_bed(
    aggregate_bed,
    aggregate_rows,
    "CIRI3_BSJ_all_samples",
    "total_bsj",
    aggregate_max,
)

sample_to_bed = dict(zip(samples, per_sample_beds))
per_sample_counts = {}
for sample in samples:
    rows = sorted(
        rows_by_sample.get(sample, []),
        key=lambda row: (
            row["chrom"],
            int(row["bsj_start"]),
            int(row["bsj_end"]),
            row["circRNA_ID"],
        ),
    )
    max_count = max((safe_float(row.get("bsj_count", 0)) for row in rows), default=0.0)
    per_sample_counts[sample] = sum(1 for row in rows if safe_float(row.get("bsj_count", 0)) > 0)
    for row in rows:
        row["sample_bsj"] = safe_float(row.get("bsj_count", 0))
    write_bed(sample_to_bed[sample], rows, f"CIRI3_BSJ_{sample}", "sample_bsj", max_count)

metadata_tsv.parent.mkdir(parents=True, exist_ok=True)
with metadata_tsv.open("w", newline="") as handle:
    fieldnames = [
        "circRNA_ID",
        "chrom",
        "bsj_start",
        "bsj_end",
        "strand",
        "gene_id",
        "circRNA_type",
        "splice_site_class",
        "abs_class",
        "total_bsj",
        "total_fsj",
        "samples_detected",
        "max_sample_bsj",
        "bsj_counts_by_sample",
        "fsj_counts_by_sample",
    ]
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    for row in aggregate_rows:
        if safe_float(row["total_bsj"]) <= 0:
            continue
        circ_id = row["circRNA_ID"]
        writer.writerow(
            {
                "circRNA_ID": circ_id,
                "chrom": row["chrom"],
                "bsj_start": row["bsj_start"],
                "bsj_end": row["bsj_end"],
                "strand": row.get("strand", "."),
                "gene_id": row.get("gene_id", "NA"),
                "circRNA_type": row.get("circRNA_type", "NA"),
                "splice_site_class": row.get("splice_site_class", "NA"),
                "abs_class": row.get("abs_class", "NA"),
                "total_bsj": format_count(row["total_bsj"]),
                "total_fsj": format_count(row["total_fsj"]),
                "samples_detected": row["samples_detected"],
                "max_sample_bsj": format_count(row["max_sample_bsj"]),
                "bsj_counts_by_sample": sample_count_string(sample_counts[circ_id], samples),
                "fsj_counts_by_sample": sample_count_string(sample_fsj_counts[circ_id], samples),
            }
        )

log_path.parent.mkdir(parents=True, exist_ok=True)
with log_path.open("w") as handle:
    aggregate_nonzero = sum(1 for row in aggregate_rows if safe_float(row["total_bsj"]) > 0)
    handle.write(f"Wrote aggregate IGV BSJ BED: {aggregate_bed}\n")
    handle.write(f"Wrote aggregate BSJ metadata TSV: {metadata_tsv}\n")
    handle.write(f"Aggregate nonzero BSJ features: {aggregate_nonzero}\n")
    for sample in samples:
        handle.write(
            f"Wrote per-sample IGV BSJ BED: {sample_to_bed[sample]} "
            f"features={per_sample_counts.get(sample, 0)}\n"
        )
