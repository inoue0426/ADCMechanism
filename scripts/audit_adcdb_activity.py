#!/usr/bin/env python3
"""Audit ADCdb cell-line activity rows by payload mechanism.

The input is a directory of ADCdb ADC detail HTML pages. The output is a flat
cell-line activity CSV plus a markdown audit report. No model is trained.
"""

from __future__ import annotations

import argparse
import csv
import html
import math
import re
from collections import Counter, defaultdict
from pathlib import Path


DETAIL_URL_BASE = "https://adcdb.idrblab.net/data/adc/details"
ACTIVITY_SECTION_TITLE = "General Information of The Activity Data Related to This ADC"
CELL_LINE_SECTION_TITLE = "Revealed Based on the Cell Line Data"
CLINICAL_SECTION_TITLE = "Identified from the Human Clinical Data"


ACTIVITY_ROW_RE = re.compile(
    r'<td class="pad-20">(?P<standard>.*?)'
    r'(?:<a href="#(?P<activity_id>[^"]+)".*?)?</td>.*?'
    r'<td class="h-center">\s*<div>\s*(?P<value>.*?)\s*</div>\s*</td>.*?'
    r'<td class="h-center">\s*<div>\s*(?P<units>.*?)\s*</div>\s*</td>.*?'
    r'<td class="pad-20">\s*<div>\s*(?P<cell_line>.*?)\s*</div>\s*</td>.*?'
    r'<td class="pad-20">\s*<div>\s*(?P<disease_model>.*?)\s*</div>\s*</td>',
    re.IGNORECASE | re.DOTALL,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--html-dir",
        type=Path,
        default=Path("data/raw/adcdb_idrblab/adc_details"),
    )
    parser.add_argument(
        "--activity-csv",
        type=Path,
        default=Path("data/processed/adcdb_cell_line_activity.csv"),
    )
    parser.add_argument(
        "--report",
        type=Path,
        default=Path("results/data_audit.md"),
    )
    return parser.parse_args()


def clean_text(text: str | None) -> str:
    if not text:
        return ""
    text = re.sub(r"<script\b.*?</script>", " ", text, flags=re.IGNORECASE | re.DOTALL)
    text = re.sub(r"<style\b.*?</style>", " ", text, flags=re.IGNORECASE | re.DOTALL)
    text = re.sub(r"<[^>]+>", " ", text)
    text = html.unescape(text).replace("\xa0", " ")
    return " ".join(text.split())


def extract_general_value(page: str, label: str) -> str:
    pattern = re.compile(
        rf"<th>\s*{re.escape(label)}\s*</th>\s*<td[^>]*>(.*?)</td>",
        re.IGNORECASE | re.DOTALL,
    )
    match = pattern.search(page)
    return clean_text(match.group(1)) if match else ""


def extract_row_link_id(page: str, label: str, href_prefix: str) -> str:
    label_pos = re.search(rf"<th>\s*{re.escape(label)}\s*</th>", page, re.IGNORECASE)
    if not label_pos:
        return ""
    window = page[label_pos.start() : label_pos.start() + 1500]
    match = re.search(re.escape(href_prefix) + r"/([^\"/]+)", window)
    return match.group(1).strip() if match else ""


def extract_cell_line_activity_block(page: str) -> str:
    activity_start = page.find(ACTIVITY_SECTION_TITLE)
    if activity_start < 0:
        return ""

    cell_line_start = page.find(CELL_LINE_SECTION_TITLE, activity_start)
    if cell_line_start < 0:
        return ""

    clinical_start = page.find(CLINICAL_SECTION_TITLE, cell_line_start)
    if clinical_start >= 0:
        return page[cell_line_start:clinical_start]

    next_unit = page.find('<div class="div-unit111">', cell_line_start + 1)
    if next_unit >= 0:
        return page[cell_line_start:next_unit]
    return page[cell_line_start:]


def split_relation_value(raw_value: str) -> tuple[str, str]:
    value = raw_value.strip()
    for relation in (">=", "<=", ">", "<", "="):
        if value.startswith(relation):
            return relation, value[len(relation) :].strip()
    return "=", value


def parse_float(value: str) -> float | None:
    match = re.search(r"[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?", value)
    if not match:
        return None
    try:
        parsed = float(match.group(0))
    except ValueError:
        return None
    return parsed if math.isfinite(parsed) else None


def is_potency_endpoint(standard_type: str) -> bool:
    text = standard_type.lower().replace(" ", "")
    return any(token in text for token in ("ic50", "lc50", "ec50", "gi50"))


def is_molar_unit(units: str) -> bool:
    return units.strip().lower() in {"pm", "nm", "um", "µm", "μm"}


def normalize_payload_identity(payload_id: str, payload_name: str) -> str:
    return payload_id or payload_name or "<missing payload>"


def mechanism_hits(payload_target: str, payload_name: str) -> set[str]:
    text = f"{payload_target} {payload_name}".lower()
    text = text.replace("topo-isomerase", "topoisomerase")

    hits: set[str] = set()
    if re.search(r"\btop(?:o|\s)?(?:isomerase)?\s*[- ]?i\b", text) or re.search(
        r"\btopoisomerase\s*1\b|\btop1\b", text
    ):
        hits.add("Topo-I inhibitor")

    microtubule_terms = [
        "microtubule",
        "tubulin",
        "auristatin",
        "monomethyl auristatin",
        "maytansinoid",
        "mertansine",
        "ansamitocin",
        "emtansine",
        "mmae",
        "mmaf",
        "dm1",
        "dm3",
        "dm4",
        "tubulysin",
        "dolastatin",
        "cryptophycin",
    ]
    if any(term in text for term in microtubule_terms):
        hits.add("microtubule-disrupting")

    dna_terms = [
        "dna",
        "calicheamicin",
        "duocarmycin",
        "duostatin",
        "pbd",
        "pyrrolobenzodiazepine",
        "tesirine",
        "indolinobenzodiazepine",
        "alkylator",
        "alkylating",
        "topoisomerase ii",
        "top2",
        "doxorubicin",
        "pnu-159682",
        "sg3199",
    ]
    topoi_terms = ["topoisomerase 1", "topoisomerase i", "top1"]
    if any(term in text for term in dna_terms) and not any(
        term in text for term in topoi_terms
    ):
        hits.add("DNA-damaging")

    return hits


def classify_mechanism(payload_target: str, payload_name: str) -> str:
    hits = mechanism_hits(payload_target, payload_name)
    if len(hits) == 1:
        return next(iter(hits))
    if len(hits) > 1:
        return "Other"
    return "Other"


def parse_page(path: Path) -> list[dict[str, str]]:
    page = path.read_text(encoding="utf-8", errors="replace")
    adc_id = path.stem
    adc_name = extract_general_value(page, "ADC Name")
    drug_status = extract_general_value(page, "Drug Status")
    antibody_name = extract_general_value(page, "Antibody Name")
    antigen_name = extract_general_value(page, "Antigen Name")
    antigen_id = extract_row_link_id(page, "Antigen Name", "/data/abt/details")
    payload_name = extract_general_value(page, "Payload Name")
    payload_id = extract_row_link_id(page, "Payload Name", "/data/payload/details")
    payload_target = extract_general_value(page, "Payload Target")
    linker_name = extract_general_value(page, "Linker Name")

    mechanism = classify_mechanism(payload_target, payload_name)
    payload_identity = normalize_payload_identity(payload_id, payload_name)
    block = extract_cell_line_activity_block(page)

    rows = []
    for row_index, match in enumerate(ACTIVITY_ROW_RE.finditer(block), start=1):
        value_raw = clean_text(match.group("value"))
        relation, value_without_relation = split_relation_value(value_raw)
        numeric_value = parse_float(value_without_relation)
        activity_id = (match.group("activity_id") or "").strip()
        rows.append(
            {
                "adc_id": adc_id,
                "adc_name": adc_name,
                "drug_status": drug_status,
                "antibody_name": antibody_name,
                "antigen_id": antigen_id,
                "antigen_name": antigen_name,
                "payload_id": payload_id,
                "payload_name": payload_name,
                "payload_identity": payload_identity,
                "payload_target": payload_target,
                "payload_mechanism": mechanism,
                "linker_name": linker_name,
                "activity_id": activity_id or f"{adc_id}_row_{row_index}",
                "standard_type": clean_text(match.group("standard")),
                "is_potency_endpoint": "1" if is_potency_endpoint(clean_text(match.group("standard"))) else "0",
                "relation": relation,
                "value_raw": value_raw,
                "value_numeric": "" if numeric_value is None else str(numeric_value),
                "units": clean_text(match.group("units")),
                "is_molar_unit": "1" if is_molar_unit(clean_text(match.group("units"))) else "0",
                "cell_line": clean_text(match.group("cell_line")),
                "disease_model": clean_text(match.group("disease_model")),
                "source_url": f"{DETAIL_URL_BASE}/{adc_id}",
                "source_html": str(path),
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "adc_id",
        "adc_name",
        "drug_status",
        "antibody_name",
        "antigen_id",
        "antigen_name",
        "payload_id",
        "payload_name",
        "payload_identity",
        "payload_target",
        "payload_mechanism",
        "linker_name",
        "activity_id",
        "standard_type",
        "is_potency_endpoint",
        "relation",
        "value_raw",
        "value_numeric",
        "units",
        "is_molar_unit",
        "cell_line",
        "disease_model",
        "source_url",
        "source_html",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def markdown_table(rows: list[dict[str, object]], headers: list[str]) -> str:
    def fmt(value: object) -> str:
        return "" if value is None else str(value)

    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(fmt(row.get(header, "")) for header in headers) + " |")
    return "\n".join(lines)


def unique_count(rows: list[dict[str, str]], column: str) -> int:
    return len({row[column] for row in rows if row[column]})


def per_mechanism_summary(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    mechanisms = sorted({row["payload_mechanism"] for row in rows})
    output = []
    for mechanism in mechanisms:
        subset = [row for row in rows if row["payload_mechanism"] == mechanism]
        unit_counter = Counter(row["units"] for row in subset)
        top_units = ", ".join(f"{unit}:{n}" for unit, n in unit_counter.most_common(4))
        payload_counter = Counter(row["payload_identity"] for row in subset)
        top_payload, top_payload_n = payload_counter.most_common(1)[0]
        output.append(
            {
                "payload_mechanism": mechanism,
                "activity_rows": len(subset),
                "nM_rows": sum(1 for row in subset if row["units"].lower() == "nm"),
                "molar_unit_rows": sum(1 for row in subset if row["is_molar_unit"] == "1"),
                "unique_payloads": unique_count(subset, "payload_identity"),
                "unique_antigens": unique_count(subset, "antigen_name"),
                "unique_cell_lines": unique_count(subset, "cell_line"),
                "unique_adcs": unique_count(subset, "adc_id"),
                "top_payload_share": f"{top_payload_n}/{len(subset)} ({top_payload_n / len(subset):.1%})",
                "top_payload": top_payload,
                "top_units": top_units,
            }
        )
    return output


def payload_summary(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    by_mechanism: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_mechanism[row["payload_mechanism"]].append(row)

    output = []
    for mechanism in sorted(by_mechanism):
        subset = by_mechanism[mechanism]
        payloads = {row["payload_identity"] for row in subset if row["payload_identity"]}
        top_payload, top_n = Counter(row["payload_identity"] for row in subset).most_common(1)[0]
        output.append(
            {
                "payload_mechanism": mechanism,
                "unique_payloads": len(payloads),
                "activity_rows": len(subset),
                "top_payload": top_payload,
                "top_payload_rows": top_n,
                "top_payload_fraction": f"{top_n / len(subset):.1%}",
            }
        )
    return output


def payload_target_summary(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    counter = Counter((row["payload_target"], row["payload_mechanism"]) for row in rows)
    output = []
    for (payload_target, mechanism), row_count in counter.most_common():
        output.append(
            {
                "payload_target": payload_target,
                "assigned_mechanism": mechanism,
                "activity_rows": row_count,
            }
        )
    return output


def cross_mechanism_cell_lines(rows: list[dict[str, str]]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    mechanisms_by_cell: dict[str, set[str]] = defaultdict(set)
    rows_by_cell: Counter[str] = Counter()
    for row in rows:
        cell_line = row["cell_line"]
        if not cell_line:
            continue
        mechanisms_by_cell[cell_line].add(row["payload_mechanism"])
        rows_by_cell[cell_line] += 1

    distribution = Counter(len(mechs) for mechs in mechanisms_by_cell.values())
    distribution_rows = [
        {"mechanism_classes_per_cell_line": key, "unique_cell_lines": distribution[key]}
        for key in sorted(distribution)
    ]

    examples = []
    for cell_line, mechanisms in mechanisms_by_cell.items():
        if len(mechanisms) <= 1:
            continue
        examples.append(
            {
                "cell_line": cell_line,
                "mechanism_count": len(mechanisms),
                "activity_rows": rows_by_cell[cell_line],
                "mechanisms": ", ".join(sorted(mechanisms)),
            }
        )
    examples.sort(key=lambda row: (-int(row["mechanism_count"]), -int(row["activity_rows"]), row["cell_line"]))
    return distribution_rows, examples


def pairwise_overlap(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    cells_by_mechanism: dict[str, set[str]] = defaultdict(set)
    for row in rows:
        if row["cell_line"]:
            cells_by_mechanism[row["payload_mechanism"]].add(row["cell_line"])
    mechanisms = sorted(cells_by_mechanism)
    output = []
    for left_idx, left in enumerate(mechanisms):
        for right in mechanisms[left_idx + 1 :]:
            output.append(
                {
                    "mechanism_a": left,
                    "mechanism_b": right,
                    "shared_cell_lines": len(cells_by_mechanism[left] & cells_by_mechanism[right]),
                }
            )
    output.sort(key=lambda row: (-int(row["shared_cell_lines"]), row["mechanism_a"], row["mechanism_b"]))
    return output


def payload_mechanism_consistency(rows: list[dict[str, str]]) -> dict[str, object]:
    mechanisms_by_payload: dict[str, set[str]] = defaultdict(set)
    rows_by_payload: Counter[str] = Counter()
    for row in rows:
        payload = row["payload_identity"]
        if not payload:
            continue
        mechanisms_by_payload[payload].add(row["payload_mechanism"])
        rows_by_payload[payload] += 1

    one_mech = sum(1 for mechanisms in mechanisms_by_payload.values() if len(mechanisms) == 1)
    total = len(mechanisms_by_payload)
    return {
        "unique_payloads": total,
        "payloads_seen_in_one_mechanism": one_mech,
        "fraction_payloads_seen_in_one_mechanism": "" if total == 0 else f"{one_mech / total:.3f}",
        "payloads_seen_in_multiple_mechanisms": total - one_mech,
    }


def write_report(
    path: Path,
    rows: list[dict[str, str]],
    html_files: list[Path],
    activity_csv: Path,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    audit_rows = [row for row in rows if row["is_potency_endpoint"] == "1"]

    per_mech = per_mechanism_summary(audit_rows)
    payloads = payload_summary(audit_rows)
    payload_targets = payload_target_summary(audit_rows)
    cell_distribution, multi_cell_examples = cross_mechanism_cell_lines(audit_rows)
    overlaps = pairwise_overlap(audit_rows)
    consistency = payload_mechanism_consistency(audit_rows)

    counts = [int(row["activity_rows"]) for row in per_mech]
    imbalance = max(counts) / min(counts) if counts and min(counts) > 0 else float("nan")

    multi_cell_count = sum(
        int(row["unique_cell_lines"])
        for row in cell_distribution
        if int(row["mechanism_classes_per_cell_line"]) > 1
    )
    total_cell_count = sum(int(row["unique_cell_lines"]) for row in cell_distribution)
    rows_in_multi_cells = sum(
        1
        for row in audit_rows
        if row["cell_line"]
        and row["cell_line"] in {example["cell_line"] for example in multi_cell_examples}
    )

    enough_cross_coverage = multi_cell_count >= 10 and rows_in_multi_cells >= 50
    no_single_payload_classes = all(int(row["unique_payloads"]) > 1 for row in payloads)
    decision = "検証可能" if enough_cross_coverage and no_single_payload_classes else "検証不能"

    causes = []
    if not enough_cross_coverage:
        causes.append("cross-mechanism cell line coverage不足")
    if not no_single_payload_classes:
        causes.append("一部mechanism classが単一payloadに近い")
    if not causes:
        causes.append("主要な即時停止理由なし")

    lines = [
        "# ADC Mechanism-Aware Response Prediction: Data Audit Pilot",
        "",
        "## Scope",
        "",
        "- Model training: not performed.",
        f"- ADCdb detail HTML pages parsed: {len(html_files):,}",
        f"- Cell-line activity rows parsed: {len(rows):,}",
        f"- Cell-line IC50/EC50/GI50 potency rows used for audit tables: {len(audit_rows):,}",
        f"- Unique ADCs with parsed potency activity: {unique_count(audit_rows, 'adc_id'):,}",
        f"- Unique payloads in potency rows: {unique_count(audit_rows, 'payload_identity'):,}",
        f"- Unique antigens in potency rows: {unique_count(audit_rows, 'antigen_name'):,}",
        f"- Unique cell lines in potency rows: {unique_count(audit_rows, 'cell_line'):,}",
        f"- Processed activity CSV: `{activity_csv}`",
        "",
        "Source scope: ADC IDs were taken from the public ADCNet `data.xlsx` file, then each ID was resolved against the public ADCdb detail page. This audits the ADCNet/ADCdb subset rather than a complete ADCdb dump.",
        "",
        "Candidate data sources checked:",
        "",
        "- Official ADCdb detail pages: chosen because they expose ADC, payload, antigen, payload target, and cell-line potency rows.",
        "- ADCNet GitHub data.xlsx: used only as the ADCdb ID seed; the table has 10 nM labels but does not retain cell-line activity rows.",
        "- AnnotationDB ADC API / digitaltumors ADCdb pipeline: useful for ADC/component metadata, but the tested individual JSON did not include cell-line activity rows.",
        "- Zenodo ADCdb mirror: identified, but direct record/file downloads returned HTTP 504 during this run.",
        "",
        "## 1. Sample Counts Per Payload Mechanism",
        "",
        markdown_table(
            per_mech,
            [
                "payload_mechanism",
                "activity_rows",
                "nM_rows",
                "molar_unit_rows",
                "unique_adcs",
                "unique_payloads",
                "unique_antigens",
                "unique_cell_lines",
                "top_payload_share",
                "top_payload",
                "top_units",
            ],
        ),
        "",
        "## 2. Unique Payloads Per Mechanism",
        "",
        markdown_table(
            payloads,
            [
                "payload_mechanism",
                "unique_payloads",
                "activity_rows",
                "top_payload",
                "top_payload_rows",
                "top_payload_fraction",
            ],
        ),
        "",
        "Payload-mechanism consistency diagnostic:",
        "",
        markdown_table([consistency], list(consistency.keys())),
        "",
        "Interpretation: payloads are expected to map mostly to one mechanism. The stronger leakage risk is when a mechanism class is represented by only one payload or is dominated by one payload.",
        "",
        "Payload target mapping used for mechanism assignment:",
        "",
        markdown_table(payload_targets, ["payload_target", "assigned_mechanism", "activity_rows"]),
        "",
        "## 3. Antigen And Cell Line Counts Per Mechanism",
        "",
        markdown_table(
            per_mech,
            ["payload_mechanism", "unique_antigens", "unique_cell_lines", "unique_adcs"],
        ),
        "",
        "## 4. Class Imbalance",
        "",
        f"- Max/min sample-count ratio across non-empty mechanism classes: {imbalance:.2f}x",
        "",
        "## 5. Cross-Mechanism Cell Line Coverage",
        "",
        f"- Cell lines measured in >1 mechanism class: {multi_cell_count:,} / {total_cell_count:,}",
        f"- Activity rows belonging to those cross-mechanism cell lines: {rows_in_multi_cells:,} / {len(audit_rows):,}",
        "",
        markdown_table(
            cell_distribution,
            ["mechanism_classes_per_cell_line", "unique_cell_lines"],
        ),
        "",
        "Pairwise shared cell-line counts:",
        "",
        markdown_table(overlaps, ["mechanism_a", "mechanism_b", "shared_cell_lines"]),
        "",
        "Top cross-mechanism cell-line examples:",
        "",
        markdown_table(
            multi_cell_examples[:15],
            ["cell_line", "mechanism_count", "activity_rows", "mechanisms"],
        ),
        "",
        "## Decision",
        "",
        f"- 判定: **{decision}**",
        f"- 主因: {', '.join(causes)}",
        "",
    ]

    if decision == "検証可能":
        lines.extend(
            [
                "Next comparison design:",
                "",
                "- Use identical train/validation/test splits for baseline and mechanism-conditioned models.",
                "- Split at least by cell line, and additionally report held-out-payload sensitivity if sample size allows.",
                "- Baseline: antigen/cell-line/delivery features without mechanism input.",
                "- Mechanism-conditioned: same features plus payload mechanism indicator or mechanism-specific head.",
                "- Primary comparison: PR-AUC/MCC/balanced accuracy on the 10 nM response label, stratified by mechanism and by held-out cell lines.",
            ]
        )
    else:
        lines.extend(
            [
                "Stop here. Do not train a model until the confounding issue is addressed by adding data or changing the audit scope.",
            ]
        )

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    html_files = sorted(args.html_dir.glob("DRG*.html"))
    rows: list[dict[str, str]] = []
    for html_path in html_files:
        rows.extend(parse_page(html_path))

    write_csv(args.activity_csv, rows)
    write_report(args.report, rows, html_files, args.activity_csv)
    print(f"Parsed {len(rows)} cell-line activity rows from {len(html_files)} HTML pages")
    print(f"Wrote {args.activity_csv}")
    print(f"Wrote {args.report}")


if __name__ == "__main__":
    main()
