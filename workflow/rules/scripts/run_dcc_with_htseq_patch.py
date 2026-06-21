#!/usr/bin/env python
"""Run DCC with a compatibility patch for newer HTSeq GTF parsing."""

from __future__ import annotations

import re
import sys


def _split_custom_exon_id(custom_exon_id: object) -> tuple[str, int, str]:
    value = str(custom_exon_id)
    transcript_id, separator, exon_text = value.rpartition(":")
    if not separator:
        raise ValueError(f"Invalid DCC custom_exon_id without exon number: {value!r}")

    match = re.match(r"(\d+)(.*)$", exon_text.strip())
    if not match:
        raise ValueError(f"Invalid DCC custom_exon_id exon number: {value!r}")

    # Newer HTSeq can keep DCC's missing final quote/semicolon in the value.
    # Preserve that suffix so adjacent-exon lookups still match the parsed keys.
    return transcript_id, int(match.group(1)), match.group(2)


def _patch_dcc_custom_exon_parser() -> None:
    from DCC.Circ_nonCirc_Exon_Match import CircNonCircExon

    def get_adjacent(self, custom_exon_id, start=True, reverse=False):
        transcript_id, exon_number, suffix = _split_custom_exon_id(custom_exon_id)
        if start:
            exon_number -= 1
        else:
            exon_number += 1
        if exon_number == -1:
            return "None"
        return f"{transcript_id}:{exon_number}{suffix}"

    def get_custom_id_num(self, custom_id):
        return _split_custom_exon_id(custom_id)[1]

    CircNonCircExon.getAdjacent = get_adjacent
    CircNonCircExon.getcustom_id_num = get_custom_id_num


def main(argv: list[str] | None = None) -> int:
    _patch_dcc_custom_exon_parser()

    from DCC.main import main as dcc_main

    sys.argv = ["DCC", *(sys.argv[1:] if argv is None else argv)]
    return dcc_main()


if __name__ == "__main__":
    raise SystemExit(main())
