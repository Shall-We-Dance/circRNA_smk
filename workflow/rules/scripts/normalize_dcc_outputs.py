import argparse
import csv
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Normalize DCC outputs into workflow BSJ matrix and annotation tables."
    )
    parser.add_argument("--samples", nargs="+", required=True)
    parser.add_argument("--circ-counts", required=True)
    parser.add_argument("--circ-coordinates", required=True)
    parser.add_argument("--linear-counts", required=True)
    parser.add_argument("--circ-skip", required=True)
    parser.add_argument("--annotation", required=True)
    parser.add_argument("--bsj", required=True)
    parser.add_argument("--fsj", required=True)
    parser.add_argument("--per-sample-dir", required=True)
    parser.add_argument("--log", required=True)
    return parser.parse_args()


def log_message(log_path: Path, message: str):
    with log_path.open("a") as handle:
        handle.write(f"[normalize_dcc_outputs] {message}\n")


def read_tsv_rows(path: Path):
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open() as handle:
        return [
            row
            for row in csv.reader(handle, delimiter="\t")
            if row and any(cell.strip() for cell in row)
        ]


def looks_like_header(row):
    if len(row) < 3:
        return True
    try:
        int(float(row[1]))
        int(float(row[2]))
    except (TypeError, ValueError):
        return True
    return False


def normalize_int(value):
    try:
        parsed = float(str(value).strip())
    except (TypeError, ValueError):
        return 0
    if parsed < 0:
        return 0
    return int(round(parsed))


def normalize_coord(value):
    try:
        parsed = int(float(str(value).strip()))
    except (TypeError, ValueError):
        return ""
    return str(parsed)


def circ_id(chrom, start, end):
    return f"{chrom}:{start}|{end}"


def locus_keys(chrom, start, end):
    chrom = str(chrom)
    keys = [(chrom, str(start), str(end))]
    if chrom.startswith("chr") and len(chrom) > 3:
        keys.append((chrom[3:], str(start), str(end)))
    else:
        keys.append((f"chr{chrom}", str(start), str(end)))
    return keys


def pick_column(header, candidates, fallback=None):
    normalized = {str(col).strip().lower(): idx for idx, col in enumerate(header)}
    for candidate in candidates:
        idx = normalized.get(candidate.lower())
        if idx is not None:
            return idx
    return fallback


def parse_circ_counts(path: Path, samples):
    rows = read_tsv_rows(path)
    if not rows:
        return []

    header = rows[0] if looks_like_header(rows[0]) else None
    data_rows = rows[1:] if header is not None else rows
    count_rows = []

    for row in data_rows:
        if len(row) < 3:
            continue
        chrom = row[0].strip()
        start = normalize_coord(row[1])
        end = normalize_coord(row[2])
        if not chrom or not start or not end:
            continue
        values = row[3:]
        counts = {
            sample: normalize_int(values[idx]) if idx < len(values) else 0
            for idx, sample in enumerate(samples)
        }
        count_rows.append(
            {
                "circRNA_ID": circ_id(chrom, start, end),
                "chrom": chrom,
                "start": start,
                "end": end,
                "counts": counts,
            }
        )
    return count_rows


def parse_coordinates(path: Path):
    rows = read_tsv_rows(path)
    if not rows:
        return {}

    default_header = [
        "chrom",
        "start",
        "end",
        "gene_id",
        "junction_type",
        "strand",
        "circ_region",
        "overall_region",
    ]
    header = rows[0] if looks_like_header(rows[0]) else default_header
    data_rows = rows[1:] if header is rows[0] else rows

    chrom_col = pick_column(header, ["chrom", "chr", "#chrom"], 0)
    start_col = pick_column(header, ["start", "chromStart"], 1)
    end_col = pick_column(header, ["end", "chromEnd"], 2)
    gene_col = pick_column(header, ["gene_id", "genename", "gene", "gene_name"], 3)
    junction_col = pick_column(header, ["junction_type", "junctiontype"], 4)
    strand_col = pick_column(header, ["strand"], 5)
    circ_region_col = pick_column(header, ["circ_region", "circRNA region"], 6)
    overall_region_col = pick_column(header, ["overall_region", "overall regions"], 7)

    annotations = {}
    for row in data_rows:
        if len(row) <= max(chrom_col, start_col, end_col):
            continue
        chrom = row[chrom_col].strip()
        start = normalize_coord(row[start_col])
        end = normalize_coord(row[end_col])
        if not chrom or not start or not end:
            continue

        def get(col, default=""):
            return row[col].strip() if col is not None and col < len(row) else default

        strand = get(strand_col, "NA")
        if strand not in {"+", "-"}:
            strand = "NA"
        record = {
            "circRNA_ID": circ_id(chrom, start, end),
            "chrom": chrom,
            "start": start,
            "end": end,
            "gene_id": get(gene_col, "NA") or "NA",
            "strand": strand,
            "junction_type": get(junction_col, "NA") or "NA",
            "circ_region": get(circ_region_col, "NA") or "NA",
            "overall_region": get(overall_region_col, "NA") or "NA",
            "source": "DCC",
        }
        for key in locus_keys(chrom, start, end):
            annotations.setdefault(key, record)
    return annotations


def ensure_raw_placeholder(path: Path, header):
    if path.exists() and path.stat().st_size > 0:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(header)


def annotation_for_count_row(row, coord_annotations):
    for key in locus_keys(row["chrom"], row["start"], row["end"]):
        if key in coord_annotations:
            record = dict(coord_annotations[key])
            record["circRNA_ID"] = row["circRNA_ID"]
            return record
    return {
        "circRNA_ID": row["circRNA_ID"],
        "chrom": row["chrom"],
        "start": row["start"],
        "end": row["end"],
        "gene_id": "NA",
        "strand": "NA",
        "junction_type": "NA",
        "circ_region": "NA",
        "overall_region": "NA",
        "source": "DCC",
    }


def write_annotation(path: Path, rows):
    fieldnames = [
        "circRNA_ID",
        "chrom",
        "start",
        "end",
        "gene_id",
        "strand",
        "junction_type",
        "circ_region",
        "overall_region",
        "source",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "NA") for field in fieldnames})


def write_matrix(path: Path, rows, samples, value_getter):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["circRNA_ID"] + list(samples))
        for row in rows:
            writer.writerow(
                [row["circRNA_ID"]]
                + [str(value_getter(row, sample)) for sample in samples]
            )


def write_per_sample_outputs(per_sample_dir: Path, samples, count_rows, annotation_rows):
    annotation_by_id = {row["circRNA_ID"]: row for row in annotation_rows}
    for sample in samples:
        sample_rows = [
            row
            for row in count_rows
            if int(row["counts"].get(sample, 0)) > 0
        ]
        sample_annotations = [
            annotation_by_id[row["circRNA_ID"]]
            for row in sample_rows
            if row["circRNA_ID"] in annotation_by_id
        ]
        prefix = per_sample_dir / f"{sample}.dcc"
        write_annotation(prefix, sample_annotations)
        write_matrix(
            prefix.with_name(prefix.name + ".BSJ_Matrix"),
            sample_rows,
            [sample],
            lambda row, sample_name: row["counts"].get(sample_name, 0),
        )
        write_matrix(
            prefix.with_name(prefix.name + ".FSJ_Matrix"),
            sample_rows,
            [sample],
            lambda row, sample_name: 0,
        )


def main():
    args = parse_args()
    samples = list(args.samples)
    circ_counts_path = Path(args.circ_counts)
    circ_coordinates_path = Path(args.circ_coordinates)
    linear_counts_path = Path(args.linear_counts)
    circ_skip_path = Path(args.circ_skip)
    log_path = Path(args.log)

    ensure_raw_placeholder(
        circ_counts_path,
        ["Chr", "Start", "End"] + samples,
    )
    ensure_raw_placeholder(
        circ_coordinates_path,
        [
            "Chr",
            "Start",
            "End",
            "genename",
            "junctiontype",
            "strand",
            "circRNA_region",
            "overall_regions",
        ],
    )
    ensure_raw_placeholder(
        linear_counts_path,
        ["Chr", "Start", "End"] + samples,
    )
    ensure_raw_placeholder(
        circ_skip_path,
        ["Chr", "Start", "End"] + samples,
    )

    count_rows = parse_circ_counts(circ_counts_path, samples)
    count_rows.sort(key=lambda row: (row["chrom"], int(row["start"]), int(row["end"])))
    coord_annotations = parse_coordinates(circ_coordinates_path)
    annotation_rows = [
        annotation_for_count_row(row, coord_annotations)
        for row in count_rows
    ]

    write_annotation(Path(args.annotation), annotation_rows)
    write_matrix(
        Path(args.bsj),
        count_rows,
        samples,
        lambda row, sample_name: row["counts"].get(sample_name, 0),
    )
    write_matrix(
        Path(args.fsj),
        count_rows,
        samples,
        lambda row, sample_name: 0,
    )
    write_per_sample_outputs(
        Path(args.per_sample_dir),
        samples,
        count_rows,
        annotation_rows,
    )
    log_message(
        log_path,
        f"Normalized {len(count_rows)} DCC circRNAs for {len(samples)} sample(s).",
    )


if __name__ == "__main__":
    main()
