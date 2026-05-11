import re
from collections import defaultdict
from pathlib import Path
from urllib.parse import quote


gtf_path = Path(snakemake.input.gtf)
gff3_path = Path(snakemake.output.gff3)
log_path = Path(snakemake.log[0])


def parse_attributes(text):
    attrs = {}
    for item in text.strip().rstrip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
        else:
            parts = item.split(None, 1)
            if len(parts) != 2:
                continue
            key, value = parts
        value = value.strip().strip('"')
        attrs[key.strip()] = value
    return attrs


def safe_id(prefix, value):
    cleaned = re.sub(r"[^A-Za-z0-9_.:-]+", "_", str(value)).strip("_")
    return f"{prefix}-{cleaned or 'unknown'}"


def fmt_attrs(attrs):
    return ";".join(
        f"{quote(str(key), safe='._:-')}={quote(str(value), safe='._:-')}"
        for key, value in attrs.items()
        if value not in (None, "")
    )


def merge_span(record, seqid, source, start, end, strand):
    if "start" not in record or "end" not in record:
        record.update(
            {
                "seqid": seqid,
                "source": source,
                "start": start,
                "end": end,
                "strand": strand,
            }
        )
    else:
        record["start"] = min(record["start"], start)
        record["end"] = max(record["end"], end)
        if record.get("strand") not in {"+", "-"} and strand in {"+", "-"}:
            record["strand"] = strand


genes = {}
transcripts = {}
exons = defaultdict(list)
rows_seen = 0

with gtf_path.open() as handle:
    for line in handle:
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) != 9:
            continue
        seqid, source, feature, start_s, end_s, score, strand, phase, attr_text = parts
        try:
            start = int(start_s)
            end = int(end_s)
        except ValueError:
            continue
        attrs = parse_attributes(attr_text)
        transcript_id = (
            attrs.get("transcript_id")
            or attrs.get("transcript")
            or attrs.get("rna_id")
            or attrs.get("ID")
        )
        gene_id = (
            attrs.get("gene_id")
            or attrs.get("gene")
            or attrs.get("gene_name")
            or attrs.get("geneID")
            or transcript_id
            or f"{seqid}:{start}-{end}"
        )
        gene_name = attrs.get("gene_name") or attrs.get("Name") or gene_id
        gene_key = gene_id
        transcript_key = transcript_id or f"{gene_id}.tx"

        gene_record = genes.setdefault(
            gene_key,
            {"id": safe_id("gene", gene_id), "name": gene_name},
        )
        transcript_record = transcripts.setdefault(
            transcript_key,
            {
                "id": safe_id("tx", transcript_key),
                "parent": gene_key,
                "name": transcript_key,
            },
        )

        if feature == "gene":
            merge_span(gene_record, seqid, source, start, end, strand)
        elif feature in {"transcript", "mRNA", "lnc_RNA", "ncRNA", "rRNA", "tRNA"}:
            merge_span(transcript_record, seqid, source, start, end, strand)
            merge_span(gene_record, seqid, source, start, end, strand)
        elif feature == "exon":
            merge_span(transcript_record, seqid, source, start, end, strand)
            merge_span(gene_record, seqid, source, start, end, strand)
            exons[transcript_key].append(
                {
                    "seqid": seqid,
                    "source": source,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "score": score,
                    "phase": phase,
                }
            )
        else:
            merge_span(gene_record, seqid, source, start, end, strand)
            merge_span(transcript_record, seqid, source, start, end, strand)
        rows_seen += 1


def sort_key(record):
    return (record.get("seqid", ""), record.get("start", 0), record.get("end", 0))


gff3_path.parent.mkdir(parents=True, exist_ok=True)
log_path.parent.mkdir(parents=True, exist_ok=True)

with gff3_path.open("w") as out:
    out.write("##gff-version 3\n")
    for gene_key, gene in sorted(genes.items(), key=lambda item: sort_key(item[1])):
        if "seqid" not in gene:
            continue
        out.write(
            "\t".join(
                [
                    gene["seqid"],
                    gene.get("source", "gtf"),
                    "gene",
                    str(gene["start"]),
                    str(gene["end"]),
                    ".",
                    gene.get("strand", "."),
                    ".",
                    fmt_attrs({"ID": gene["id"], "Name": gene.get("name", gene_key)}),
                ]
            )
            + "\n"
        )
        gene_transcripts = [
            (tx_key, tx)
            for tx_key, tx in transcripts.items()
            if tx.get("parent") == gene_key and "seqid" in tx
        ]
        for tx_key, tx in sorted(gene_transcripts, key=lambda item: sort_key(item[1])):
            out.write(
                "\t".join(
                    [
                        tx["seqid"],
                        tx.get("source", "gtf"),
                        "mRNA",
                        str(tx["start"]),
                        str(tx["end"]),
                        ".",
                        tx.get("strand", "."),
                        ".",
                        fmt_attrs(
                            {
                                "ID": tx["id"],
                                "Parent": gene["id"],
                                "Name": tx.get("name", tx_key),
                            }
                        ),
                    ]
                )
                + "\n"
            )
            for index, exon in enumerate(
                sorted(exons.get(tx_key, []), key=lambda rec: (rec["start"], rec["end"])),
                start=1,
            ):
                out.write(
                    "\t".join(
                        [
                            exon["seqid"],
                            exon.get("source", "gtf"),
                            "exon",
                            str(exon["start"]),
                            str(exon["end"]),
                            ".",
                            exon.get("strand", "."),
                            ".",
                            fmt_attrs(
                                {
                                    "ID": f"{tx['id']}.exon{index}",
                                    "Parent": tx["id"],
                                }
                            ),
                        ]
                    )
                    + "\n"
                )

log_path.write_text(
    f"Converted {rows_seen} GTF rows to {gff3_path} for rmats2sashimiplot coordinate mode.\n"
)
