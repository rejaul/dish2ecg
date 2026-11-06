"""Recover tabular text from Enhanced Metafile (EMF) images.

Some supplementary tables are published as vector images pasted into a Word
document rather than as a table or a spreadsheet. Supplemental Table 1 of the
CiPA MEA pilot study (Millard D et al., *Toxicol Sci* 2018;164:550-562,
doi:10.1093/toxsci/kfy110) is one: the mean baseline beat period, field
potential duration and rate-corrected field potential duration of all 18
participating studies exist only inside `word/media/image1.emf`.

EMF stores drawn text in EMR_EXTTEXTOUTW (record type 84) and EMR_EXTTEXTOUTA
(83) records together with its reference point, so the strings and their
positions are recoverable exactly - this is reading the file format, not
optical character recognition, and the digits come back verbatim.
"""
from __future__ import annotations

import collections
import struct
import zipfile
from pathlib import Path

__all__ = ["emf_text_records", "emf_table", "docx_emf_parts"]

_EMR_EXTTEXTOUTA = 83
_EMR_EXTTEXTOUTW = 84


def emf_text_records(data: bytes) -> list[tuple[int, int, str]]:
    """Every text record in an EMF blob as ``(y, x, string)``, document order."""
    out: list[tuple[int, int, str]] = []
    i, n = 0, len(data)
    while i + 8 <= n:
        rtype, rsize = struct.unpack_from("<II", data, i)
        if rsize < 8 or i + rsize > n:
            break
        if rtype in (_EMR_EXTTEXTOUTA, _EMR_EXTTEXTOUTW):
            try:
                # record header (8) + rclBounds (16) + iGraphicsMode (4)
                # + exScale (4) + eyScale (4), then the EmrText structure
                off = i + 8 + 16 + 4 + 4 + 4
                ref_x, ref_y = struct.unpack_from("<ii", data, off)
                off += 8
                chars = struct.unpack_from("<I", data, off)[0]
                off += 4
                off_string = struct.unpack_from("<I", data, off)[0]
                start = i + off_string
                width = 2 if rtype == _EMR_EXTTEXTOUTW else 1
                if 0 < chars < 4000 and start + chars * width <= n:
                    raw = data[start:start + chars * width]
                    s = (raw.decode("utf-16-le", "replace") if width == 2
                         else raw.decode("latin-1", "replace"))
                    out.append((ref_y, ref_x, s))
            except struct.error:
                pass
        i += rsize
    return out


def emf_table(data: bytes, *, row_tolerance: int = 6) -> list[list[str]]:
    """Reconstruct a table by grouping text records into rows by y, then x.

    ``row_tolerance`` merges reference points within that many units into one
    row, for tables whose cells are not drawn on a perfectly shared baseline.
    """
    recs = emf_text_records(data)
    if not recs:
        return []
    buckets: dict[int, list[tuple[int, str]]] = collections.defaultdict(list)
    anchors: list[int] = []
    for y, x, s in recs:
        hit = next((a for a in anchors if abs(a - y) <= row_tolerance), None)
        if hit is None:
            anchors.append(y)
            hit = y
        buckets[hit].append((x, s))
    return [[s for _, s in sorted(buckets[y])] for y in sorted(buckets)]


def docx_emf_parts(path: str | Path) -> dict[str, bytes]:
    """Every ``word/media/*.emf`` part of a .docx, keyed by part name."""
    with zipfile.ZipFile(path) as z:
        return {n: z.read(n) for n in z.namelist()
                if n.lower().startswith("word/media/") and n.lower().endswith(".emf")}
