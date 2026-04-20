from pathlib import Path
import csv

pairwise_files = [Path(p) for p in snakemake.input.results]
pairs = list(snakemake.params.pairs)
out_path = Path(snakemake.output.all_samples)

if len(pairwise_files) != len(pairs):
    raise ValueError("Number of pairwise result files does not match number of pair names")

header = None
rows_out = []
for pair, file_path in zip(pairs, pairwise_files):
    with file_path.open() as fh:
        reader = csv.reader(fh, delimiter="\t")
        file_header = next(reader, None)
        if not file_header:
            continue
        if header is None:
            header = file_header
        elif file_header != header:
            raise ValueError(
                f"Header mismatch while merging '{file_path}'.\n"
                f"Expected: {header}\nGot: {file_header}"
            )

        for row in reader:
            if not row:
                continue
            rows_out.append([pair] + row)

out_path.parent.mkdir(parents=True, exist_ok=True)
with out_path.open("w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t")
    if header is None:
        writer.writerow(["comparison"])
    else:
        writer.writerow(["comparison"] + header)
        writer.writerows(rows_out)
