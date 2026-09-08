#!/usr/bin/env python3
"""Fetch ADCdb detail pages for ADC IDs listed in an xlsx file.

This script only downloads public HTML detail pages. It does not train models or
derive labels.
"""

from __future__ import annotations

import argparse
import csv
import html
import re
import time
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path
from zipfile import ZipFile


BASE_URL = "https://adcdb.idrblab.net/data/adc/details"
XLSX_NS = {"a": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--xlsx",
        type=Path,
        default=Path("data/raw/adcnet_github/data.xlsx"),
        help="xlsx file containing an 'ADC ID' column",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("data/raw/adcdb_idrblab/adc_details"),
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("data/raw/adcdb_idrblab/adc_details_manifest.csv"),
    )
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--delay-sec", type=float, default=0.25)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def column_index(cell_ref: str) -> int:
    match = re.match(r"[A-Z]+", cell_ref)
    if not match:
        raise ValueError(f"Unexpected cell reference: {cell_ref}")
    value = 0
    for char in match.group(0):
        value = value * 26 + ord(char) - ord("A") + 1
    return value - 1


def read_xlsx_rows(path: Path) -> list[list[str]]:
    with ZipFile(path) as archive:
        shared_strings: list[str] = []
        if "xl/sharedStrings.xml" in archive.namelist():
            shared_root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            for item in shared_root.findall("a:si", XLSX_NS):
                shared_strings.append(
                    "".join(text.text or "" for text in item.findall(".//a:t", XLSX_NS))
                )

        sheet_root = ET.fromstring(archive.read("xl/worksheets/sheet1.xml"))
        rows: list[list[str]] = []
        for row in sheet_root.findall(".//a:sheetData/a:row", XLSX_NS):
            values: list[str] = []
            for cell in row.findall("a:c", XLSX_NS):
                idx = column_index(cell.attrib["r"])
                while len(values) <= idx:
                    values.append("")

                cell_type = cell.attrib.get("t")
                value_node = cell.find("a:v", XLSX_NS)
                if cell_type == "s" and value_node is not None:
                    values[idx] = shared_strings[int(value_node.text or "0")]
                elif cell_type == "inlineStr":
                    values[idx] = "".join(
                        text.text or "" for text in cell.findall(".//a:t", XLSX_NS)
                    )
                elif value_node is not None:
                    values[idx] = html.unescape(value_node.text or "")
            rows.append(values)
        return rows


def read_adc_ids(path: Path) -> list[str]:
    rows = read_xlsx_rows(path)
    if not rows:
        return []
    header = rows[0]
    try:
        adc_id_idx = header.index("ADC ID")
    except ValueError as exc:
        raise ValueError(f"Missing 'ADC ID' column in {path}") from exc

    ids = []
    seen = set()
    for row in rows[1:]:
        adc_id = row[adc_id_idx].strip() if adc_id_idx < len(row) else ""
        if adc_id and adc_id not in seen:
            ids.append(adc_id)
            seen.add(adc_id)
    return ids


def fetch_page(adc_id: str, out_path: Path) -> tuple[str, int, str]:
    url = f"{BASE_URL}/{adc_id}"
    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": "ADCMechanism-data-audit/0.1 (+https://github.com/inoue0426/ADCMechanism)"
        },
    )
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            body = response.read()
            out_path.write_bytes(body)
            return "ok", len(body), ""
    except urllib.error.HTTPError as exc:
        return f"http_{exc.code}", 0, str(exc)
    except urllib.error.URLError as exc:
        return "url_error", 0, str(exc.reason)
    except TimeoutError as exc:
        return "timeout", 0, str(exc)


def main() -> None:
    args = parse_args()
    adc_ids = read_adc_ids(args.xlsx)
    if args.limit is not None:
        adc_ids = adc_ids[: args.limit]

    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.manifest.parent.mkdir(parents=True, exist_ok=True)

    records = []
    for index, adc_id in enumerate(adc_ids, start=1):
        out_path = args.out_dir / f"{adc_id}.html"
        if out_path.exists() and out_path.stat().st_size > 0 and not args.force:
            status = "cached"
            size = out_path.stat().st_size
            error = ""
        else:
            status, size, error = fetch_page(adc_id, out_path)
            time.sleep(args.delay_sec)
        records.append(
            {
                "adc_id": adc_id,
                "url": f"{BASE_URL}/{adc_id}",
                "path": str(out_path),
                "status": status,
                "bytes": size,
                "error": error,
            }
        )
        if index % 25 == 0 or index == len(adc_ids):
            print(f"{index}/{len(adc_ids)} pages processed")

    with args.manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["adc_id", "url", "path", "status", "bytes", "error"]
        )
        writer.writeheader()
        writer.writerows(records)

    ok = sum(row["status"] in {"ok", "cached"} for row in records)
    print(f"Wrote {args.manifest} with {ok}/{len(records)} available pages")


if __name__ == "__main__":
    main()
