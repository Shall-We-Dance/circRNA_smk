#!/usr/bin/env python3
"""Write explicit empty outputs for circtools circtest no-candidate runs."""

import argparse
import csv
import textwrap
from pathlib import Path
from xml.sax.saxutils import escape
from zipfile import ZIP_DEFLATED, ZipFile


NO_CANDIDATES_MESSAGE = (
    "circtools circtest completed filtering/statistical testing, but no "
    "candidates passed the configured thresholds, so circtools did not emit "
    "its normal result files."
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create explicit empty circtest CSV/XLSX/PDF outputs."
    )
    parser.add_argument("--csv", required=True, help="Output circtest.csv path.")
    parser.add_argument("--xlsx", required=True, help="Output circtest.xlsx path.")
    parser.add_argument("--pdf", required=True, help="Output circtest.pdf path.")
    parser.add_argument("--comparison", required=True, help="Comparison name.")
    parser.add_argument("--log", required=True, help="circtest log path.")
    return parser.parse_args()


def ensure_parent(path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)


def write_csv(path, comparison, log_path):
    ensure_parent(path)
    with Path(path).open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["status", "comparison", "message", "source_log"])
        writer.writerow(["no_candidates", comparison, NO_CANDIDATES_MESSAGE, log_path])


def column_name(index):
    name = ""
    index += 1
    while index:
        index, remainder = divmod(index - 1, 26)
        name = chr(ord("A") + remainder) + name
    return name


def worksheet_xml(rows):
    row_xml = []
    for row_index, row in enumerate(rows, start=1):
        cells = []
        for col_index, value in enumerate(row):
            ref = f"{column_name(col_index)}{row_index}"
            text = escape(str(value))
            cells.append(
                f'<c r="{ref}" t="inlineStr"><is><t>{text}</t></is></c>'
            )
        row_xml.append(f'<row r="{row_index}">{"".join(cells)}</row>')
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        f'<sheetData>{"".join(row_xml)}</sheetData>'
        "</worksheet>"
    )


def write_xlsx(path, comparison, log_path):
    ensure_parent(path)
    summary_rows = [
        ["status", "comparison", "message", "source_log"],
        ["no_candidates", comparison, NO_CANDIDATES_MESSAGE, log_path],
    ]
    count_rows = [
        ["status", "message"],
        ["not_applicable", "No significant circles were available for this sheet."],
    ]
    sheets = [
        ("Significant circles", summary_rows),
        ("Circle Counts", count_rows),
        ("Linear Counts", count_rows),
    ]

    content_types = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
  <Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
  <Default Extension="xml" ContentType="application/xml"/>
  <Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>
  <Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>
  <Override PartName="/xl/worksheets/sheet2.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>
  <Override PartName="/xl/worksheets/sheet3.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>
  <Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.core-properties+xml"/>
  <Override PartName="/docProps/app.xml" ContentType="application/vnd.openxmlformats-officedocument.extended-properties+xml"/>
</Types>
"""
    rels = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>
  <Relationship Id="rId2" Type="http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties" Target="docProps/core.xml"/>
  <Relationship Id="rId3" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/extended-properties" Target="docProps/app.xml"/>
</Relationships>
"""
    workbook_sheets = "".join(
        f'<sheet name="{escape(name)}" sheetId="{i}" r:id="rId{i}"/>'
        for i, (name, _) in enumerate(sheets, start=1)
    )
    workbook = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
        f"<sheets>{workbook_sheets}</sheets>"
        "</workbook>"
    )
    workbook_rels = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/>
  <Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet2.xml"/>
  <Relationship Id="rId3" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet3.xml"/>
</Relationships>
"""
    core = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" xmlns:dc="http://purl.org/dc/elements/1.1/">
  <dc:title>Empty circtools circtest result</dc:title>
</cp:coreProperties>
"""
    app = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties">
  <Application>circRNA_smk</Application>
</Properties>
"""

    with ZipFile(path, "w", ZIP_DEFLATED) as zf:
        zf.writestr("[Content_Types].xml", content_types)
        zf.writestr("_rels/.rels", rels)
        zf.writestr("xl/workbook.xml", workbook)
        zf.writestr("xl/_rels/workbook.xml.rels", workbook_rels)
        zf.writestr("docProps/core.xml", core)
        zf.writestr("docProps/app.xml", app)
        for i, (_, rows) in enumerate(sheets, start=1):
            zf.writestr(f"xl/worksheets/sheet{i}.xml", worksheet_xml(rows))


def pdf_escape(text):
    return str(text).replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")


def write_pdf(path, comparison, log_path):
    ensure_parent(path)
    raw_lines = [
        "circtools circtest: no candidates",
        f"Comparison: {comparison}",
        NO_CANDIDATES_MESSAGE,
        f"Source log: {log_path}",
    ]
    lines = [raw_lines[0]]
    for line in raw_lines[1:]:
        lines.extend(textwrap.wrap(line, width=86) or [""])

    text_ops = ["BT", "/F1 16 Tf", "72 720 Td", f"({pdf_escape(lines[0])}) Tj"]
    text_ops.extend(["/F1 10 Tf"])
    for line in lines[1:]:
        text_ops.append("0 -24 Td")
        text_ops.append(f"({pdf_escape(line)}) Tj")
    text_ops.append("ET")
    stream = "\n".join(text_ops).encode("utf-8")

    objects = [
        b"<< /Type /Catalog /Pages 2 0 R >>",
        b"<< /Type /Pages /Kids [3 0 R] /Count 1 >>",
        (
            b"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 612 792] "
            b"/Resources << /Font << /F1 4 0 R >> >> /Contents 5 0 R >>"
        ),
        b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>",
        b"<< /Length " + str(len(stream)).encode("ascii") + b" >>\nstream\n"
        + stream + b"\nendstream",
    ]

    chunks = [b"%PDF-1.4\n"]
    offsets = []
    for obj_id, obj in enumerate(objects, start=1):
        offsets.append(sum(len(chunk) for chunk in chunks))
        chunks.append(f"{obj_id} 0 obj\n".encode("ascii"))
        chunks.append(obj)
        chunks.append(b"\nendobj\n")
    xref_offset = sum(len(chunk) for chunk in chunks)
    chunks.append(f"xref\n0 {len(objects) + 1}\n".encode("ascii"))
    chunks.append(b"0000000000 65535 f \n")
    for offset in offsets:
        chunks.append(f"{offset:010d} 00000 n \n".encode("ascii"))
    chunks.append(
        (
            f"trailer\n<< /Size {len(objects) + 1} /Root 1 0 R >>\n"
            f"startxref\n{xref_offset}\n%%EOF\n"
        ).encode("ascii")
    )
    Path(path).write_bytes(b"".join(chunks))


def append_log(log_path, outputs):
    with Path(log_path).open("a") as handle:
        handle.write("\n[circRNA_smk] No circtest candidates detected; ")
        handle.write("wrote explicit empty outputs:\n")
        for output in outputs:
            handle.write(f"[circRNA_smk]   {output}\n")


def main():
    args = parse_args()
    write_csv(args.csv, args.comparison, args.log)
    write_xlsx(args.xlsx, args.comparison, args.log)
    write_pdf(args.pdf, args.comparison, args.log)
    append_log(args.log, [args.csv, args.xlsx, args.pdf])


if __name__ == "__main__":
    main()
