#!/usr/bin/env python3
"""ADC payload-mechanism identifiability audit and cheap falsification.

This script intentionally uses only standard-library code. It reads the local
raw/processed ADCdb artifacts, derives a conservative 10 nM cell-line response
dataset, audits mechanism confounding, then runs small categorical baselines.

No neural model is trained.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
import statistics
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from zipfile import ZipFile


XLSX_NS = {
    "main": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
    "rel": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
    "pkgrel": "http://schemas.openxmlformats.org/package/2006/relationships",
}

MAJOR_MECHANISMS = {
    "DNA-damaging",
    "Topo-I inhibitor",
    "microtubule-disrupting",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-dir", type=Path, default=Path("data/raw"))
    parser.add_argument(
        "--activity-csv",
        type=Path,
        default=Path("data/processed/adcdb_cell_line_activity.csv"),
    )
    parser.add_argument(
        "--context-csv",
        type=Path,
        default=Path("data/processed/adc_response_10nM_context.csv"),
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("results"),
    )
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--seed", type=int, default=44)
    parser.add_argument("--shuffle-repeats", type=int, default=50)
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        keys: list[str] = []
        seen: set[str] = set()
        for row in rows:
            for key in row:
                if key not in seen:
                    seen.add(key)
                    keys.append(key)
        fieldnames = keys
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def markdown_table(rows: list[dict[str, object]], headers: list[str]) -> str:
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(str(row.get(header, "")) for header in headers) + " |")
    return "\n".join(lines)


def col_to_idx(cell_ref: str) -> int:
    letters = "".join(char for char in cell_ref if char.isalpha())
    out = 0
    for char in letters:
        out = out * 26 + ord(char.upper()) - 64
    return out - 1


def load_shared_strings(zipped: ZipFile) -> list[str]:
    try:
        root = ET.fromstring(zipped.read("xl/sharedStrings.xml"))
    except KeyError:
        return []
    shared: list[str] = []
    for si in root.findall("main:si", XLSX_NS):
        parts = [
            node.text or ""
            for node in si.iter("{http://schemas.openxmlformats.org/spreadsheetml/2006/main}t")
        ]
        shared.append("".join(parts))
    return shared


def cell_value(cell: ET.Element, shared: list[str]) -> str:
    cell_type = cell.attrib.get("t")
    if cell_type == "inlineStr":
        return "".join(
            node.text or ""
            for node in cell.iter("{http://schemas.openxmlformats.org/spreadsheetml/2006/main}t")
        )

    value = cell.find("main:v", XLSX_NS)
    if value is None:
        return ""
    raw = value.text or ""
    if cell_type == "s":
        try:
            return shared[int(raw)]
        except (ValueError, IndexError):
            return raw
    return raw


def read_xlsx_first_sheet(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with ZipFile(path) as zipped:
        shared = load_shared_strings(zipped)
        root = ET.fromstring(zipped.read("xl/worksheets/sheet1.xml"))
        rows: list[list[str]] = []
        for row in root.findall(".//main:sheetData/main:row", XLSX_NS):
            values: list[str] = []
            for cell in row.findall("main:c", XLSX_NS):
                index = col_to_idx(cell.attrib.get("r", "A1"))
                while len(values) <= index:
                    values.append("")
                values[index] = cell_value(cell, shared)
            if any(value.strip() for value in values):
                rows.append(values)

    if not rows:
        return [], []
    headers = rows[0]
    records: list[dict[str, str]] = []
    for row in rows[1:]:
        padded = row + [""] * max(0, len(headers) - len(row))
        records.append(dict(zip(headers, padded[: len(headers)])))
    return headers, records


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def summarize_raw_files(raw_dir: Path) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    files = sorted(path for path in raw_dir.rglob("*") if path.is_file())
    inventory: list[dict[str, object]] = []
    for path in files:
        rel = path.relative_to(raw_dir)
        top = rel.parts[0] if rel.parts else "."
        note = ""
        if rel.as_posix() == "adcnet_github/data.xlsx":
            note = "ADCNet ADC-level table with ADC IDs, sequences, SMILES, DAR, and threshold labels."
        elif rel.as_posix() == "adcnet_github/files_data.xlsx":
            note = "Byte-identical copy of adcnet_github/data.xlsx."
        elif rel.as_posix() == "adcnet_github/t_data.xlsx":
            note = "ADCNet auxiliary table; has fewer rows and lacks a 10 nM label column."
        elif rel.as_posix() == "adcnet_github/files_api.json":
            note = "GitHub API listing for ADCNet files."
        elif rel.as_posix() == "adcdb_idrblab/adc_details_manifest.csv":
            note = "Manifest for locally cached ADCdb detail pages."
        elif rel.as_posix().startswith("adcdb_idrblab/adc_details/"):
            note = "Official ADCdb ADC detail HTML with component metadata and cell-line activity rows."
        elif rel.as_posix().startswith("digitaltumors_adcdb_pipeline/"):
            note = "External pipeline documentation/config; useful for metadata provenance, not response rows here."
        elif rel.as_posix().startswith("annotationdb_adcdb/"):
            note = "AnnotationDB ADC metadata snapshot/API schema; local sample does not expose cell-line response rows."
        elif rel.suffix == ".html":
            note = "ADCdb HTML page cached during source inspection."
        inventory.append(
            {
                "path": str(path),
                "relative_path": rel.as_posix(),
                "top_level": top,
                "suffix": rel.suffix.lower() or "<none>",
                "bytes": path.stat().st_size,
                "sha256": file_sha256(path),
                "note": note,
            }
        )

    summary: list[dict[str, object]] = []
    by_top: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in inventory:
        by_top[str(row["top_level"])].append(row)
    for top, rows in sorted(by_top.items()):
        summary.append(
            {
                "top_level": top,
                "files": len(rows),
                "total_bytes": sum(int(row["bytes"]) for row in rows),
                "suffixes": ", ".join(
                    f"{suffix}:{count}"
                    for suffix, count in sorted(Counter(str(row["suffix"]) for row in rows).items())
                ),
            }
        )
    return inventory, summary


def summarize_adcnet_xlsx(raw_dir: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for path in sorted((raw_dir / "adcnet_github").glob("*.xlsx")):
        headers, records = read_xlsx_first_sheet(path)
        label_columns = [header for header in headers if "label" in header.lower()]
        label_summaries = []
        for column in label_columns:
            counts = Counter(record.get(column, "") for record in records)
            label_summaries.append(
                column
                + ":"
                + ",".join(f"{key or '<blank>'}={counts[key]}" for key in sorted(counts))
            )
        rows.append(
            {
                "path": str(path),
                "rows": len(records),
                "columns": len(headers),
                "headers": "; ".join(headers),
                "label_columns": "; ".join(label_columns),
                "label_counts": "; ".join(label_summaries),
                "sha256": file_sha256(path),
            }
        )
    return rows


def to_nm(value: str, units: str) -> float | None:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return None
    unit = units.strip().lower()
    if unit == "pm":
        return numeric / 1000.0
    if unit == "nm":
        return numeric
    if unit in {"um", "umol/l", "µm", "μm"}:
        return numeric * 1000.0
    return None


def deterministic_label(relation: str, value_nm: float | None, threshold_nm: float = 10.0) -> int | None:
    if value_nm is None or not math.isfinite(value_nm):
        return None
    if relation == "=":
        return 1 if value_nm <= threshold_nm else 0
    if relation in {"<", "<="}:
        return 1 if value_nm <= threshold_nm else None
    if relation == ">":
        return 0 if value_nm >= threshold_nm else None
    if relation == ">=":
        return 0 if value_nm > threshold_nm else None
    return None


def is_potency(row: dict[str, str]) -> bool:
    if row.get("is_potency_endpoint") == "1":
        return True
    text = row.get("standard_type", "").lower().replace(" ", "")
    return any(token in text for token in ("ic50", "lc50", "ec50", "gi50"))


def derive_usable_rows(rows: list[dict[str, str]]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    usable: list[dict[str, object]] = []
    exclusions: list[dict[str, object]] = []
    for row in rows:
        if not is_potency(row):
            reason = "non_potency_endpoint"
        else:
            value_nm = to_nm(row.get("value_numeric", ""), row.get("units", ""))
            if value_nm is None:
                reason = "non_molar_or_missing_value"
            else:
                label = deterministic_label(row.get("relation", ""), value_nm)
                if label is None:
                    reason = "ambiguous_censored_value"
                else:
                    out = dict(row)
                    out["value_nM"] = f"{value_nm:.12g}"
                    out["label_10nM"] = str(label)
                    usable.append(out)
                    continue

        exclusions.append(
            {
                "adc_id": row.get("adc_id", ""),
                "activity_id": row.get("activity_id", ""),
                "standard_type": row.get("standard_type", ""),
                "relation": row.get("relation", ""),
                "value_numeric": row.get("value_numeric", ""),
                "units": row.get("units", ""),
                "payload_mechanism": row.get("payload_mechanism", ""),
                "reason": reason,
            }
        )
    return usable, exclusions


def collapse_context_rows(rows: list[dict[str, object]]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    grouped: dict[tuple[str, str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[
            (
                str(row["adc_id"]),
                str(row["cell_line"]),
                str(row["antigen_name"]),
                str(row["payload_identity"]),
            )
        ].append(row)

    collapsed: list[dict[str, object]] = []
    conflicts: list[dict[str, object]] = []
    for group_rows in grouped.values():
        best = min(group_rows, key=lambda row: float(row["value_nM"]))
        labels = {str(row["label_10nM"]) for row in group_rows}
        value_min = min(float(row["value_nM"]) for row in group_rows)
        output = dict(best)
        output["value_nM_min"] = f"{value_min:.12g}"
        output["label_10nM"] = "1" if "1" in labels else "0"
        output["source_activity_rows"] = str(len(group_rows))
        output["has_conflicting_10nM_labels"] = "1" if len(labels) > 1 else "0"
        output["source_activity_ids"] = ";".join(str(row["activity_id"]) for row in group_rows)
        collapsed.append(output)

        if len(labels) > 1:
            conflicts.append(
                {
                    "adc_id": best["adc_id"],
                    "cell_line": best["cell_line"],
                    "antigen_name": best["antigen_name"],
                    "payload_identity": best["payload_identity"],
                    "payload_mechanism": best["payload_mechanism"],
                    "labels": ",".join(sorted(labels)),
                    "source_activity_rows": len(group_rows),
                    "values_nM": ";".join(str(row["value_nM"]) for row in group_rows),
                    "activity_ids": ";".join(str(row["activity_id"]) for row in group_rows),
                }
            )

    collapsed.sort(
        key=lambda row: (
            str(row["payload_mechanism"]),
            str(row["payload_identity"]),
            str(row["antigen_name"]),
            str(row["cell_line"]),
            str(row["adc_id"]),
        )
    )
    return collapsed, conflicts


def unique_count(rows: list[dict[str, object]], key: str) -> int:
    return len({str(row.get(key, "")) for row in rows if str(row.get(key, ""))})


def prevalence(rows: list[dict[str, object]]) -> float:
    if not rows:
        return float("nan")
    return sum(int(row["label_10nM"]) for row in rows) / len(rows)


def entropy_effective_n(counter: Counter[str]) -> float:
    total = sum(counter.values())
    if total == 0:
        return 0.0
    entropy = 0.0
    for count in counter.values():
        p = count / total
        entropy -= p * math.log(p)
    return math.exp(entropy)


def per_mechanism_table(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    for mechanism in sorted({str(row["payload_mechanism"]) for row in rows}):
        subset = [row for row in rows if row["payload_mechanism"] == mechanism]
        payload_counts = Counter(str(row["payload_identity"]) for row in subset)
        antigen_counts = Counter(str(row["antigen_name"]) for row in subset)
        cell_counts = Counter(str(row["cell_line"]) for row in subset)
        top_payload, top_payload_n = payload_counts.most_common(1)[0]
        top_antigen, top_antigen_n = antigen_counts.most_common(1)[0]
        top_cell, top_cell_n = cell_counts.most_common(1)[0]
        output.append(
            {
                "payload_mechanism": mechanism,
                "observations": len(subset),
                "positives_10nM": sum(int(row["label_10nM"]) for row in subset),
                "prevalence_10nM": f"{prevalence(subset):.3f}",
                "unique_payloads": len(payload_counts),
                "effective_payload_n": f"{entropy_effective_n(payload_counts):.2f}",
                "unique_antigens": len(antigen_counts),
                "unique_cell_lines": len(cell_counts),
                "top_payload": top_payload,
                "top_payload_share": f"{top_payload_n / len(subset):.1%}",
                "top_antigen": top_antigen,
                "top_antigen_share": f"{top_antigen_n / len(subset):.1%}",
                "top_cell_line": top_cell,
                "top_cell_share": f"{top_cell_n / len(subset):.1%}",
            }
        )
    return output


def matrix_counts(rows: list[dict[str, object]], row_key: str, col_key: str) -> list[dict[str, object]]:
    row_values = sorted({str(row[row_key]) for row in rows})
    col_values = sorted({str(row[col_key]) for row in rows})
    counts = Counter((str(row[row_key]), str(row[col_key])) for row in rows)
    output = []
    for row_value in row_values:
        line: dict[str, object] = {row_key: row_value}
        for col_value in col_values:
            line[col_value] = counts[(row_value, col_value)]
        output.append(line)
    return output


def cramers_v(rows: list[dict[str, object]], key_a: str, key_b: str) -> float:
    values_a = sorted({str(row[key_a]) for row in rows})
    values_b = sorted({str(row[key_b]) for row in rows})
    if len(values_a) < 2 or len(values_b) < 2:
        return float("nan")
    counts = Counter((str(row[key_a]), str(row[key_b])) for row in rows)
    row_totals = Counter(str(row[key_a]) for row in rows)
    col_totals = Counter(str(row[key_b]) for row in rows)
    total = len(rows)
    chi2 = 0.0
    for value_a in values_a:
        for value_b in values_b:
            expected = row_totals[value_a] * col_totals[value_b] / total
            if expected > 0:
                observed = counts[(value_a, value_b)]
                chi2 += (observed - expected) ** 2 / expected
    phi2 = chi2 / total
    denom = min(len(values_a) - 1, len(values_b) - 1)
    return math.sqrt(phi2 / denom) if denom > 0 else float("nan")


def association_table(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    definitions = [
        ("payload_identity", "mechanism_payload"),
        ("antigen_name", "mechanism_antigen"),
        ("cell_line", "mechanism_cell_line"),
        ("disease_model", "mechanism_disease"),
        ("cell_antigen_context", "mechanism_cell_antigen_context"),
    ]
    output = []
    enriched = [dict(row, cell_antigen_context=f"{row['cell_line']} || {row['antigen_name']}") for row in rows]
    for key, label in definitions:
        output.append(
            {
                "association": label,
                "left": "payload_mechanism",
                "right": key,
                "cramers_v": f"{cramers_v(enriched, 'payload_mechanism', key):.3f}",
                "right_unique_values": unique_count(enriched, key),
            }
        )
    return output


def overlap_table(rows: list[dict[str, object]], context_key: str) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    mechanisms_by_context: dict[str, set[str]] = defaultdict(set)
    payloads_by_context: dict[str, set[str]] = defaultdict(set)
    labels_by_context: dict[str, set[str]] = defaultdict(set)
    rows_by_context: Counter[str] = Counter()
    for row in rows:
        context = str(row[context_key])
        mechanisms_by_context[context].add(str(row["payload_mechanism"]))
        payloads_by_context[context].add(str(row["payload_identity"]))
        labels_by_context[context].add(str(row["label_10nM"]))
        rows_by_context[context] += 1

    distribution = Counter(len(mechanisms) for mechanisms in mechanisms_by_context.values())
    summary = [
        {
            "context": context_key,
            "mechanism_count": degree,
            "unique_contexts": distribution[degree],
        }
        for degree in sorted(distribution)
    ]

    examples: list[dict[str, object]] = []
    for context, mechanisms in mechanisms_by_context.items():
        if len(mechanisms) <= 1:
            continue
        examples.append(
            {
                "context": context,
                "mechanism_count": len(mechanisms),
                "observations": rows_by_context[context],
                "unique_payloads": len(payloads_by_context[context]),
                "has_label_variation": "1" if len(labels_by_context[context]) > 1 else "0",
                "mechanisms": "; ".join(sorted(mechanisms)),
            }
        )
    examples.sort(
        key=lambda row: (
            -int(row["mechanism_count"]),
            -int(row["observations"]),
            str(row["context"]),
        )
    )
    return summary, examples


def pairwise_context_overlap(rows: list[dict[str, object]], context_key: str) -> list[dict[str, object]]:
    contexts_by_mechanism: dict[str, set[str]] = defaultdict(set)
    for row in rows:
        contexts_by_mechanism[str(row["payload_mechanism"])].add(str(row[context_key]))
    mechanisms = sorted(contexts_by_mechanism)
    output: list[dict[str, object]] = []
    for left_index, left in enumerate(mechanisms):
        for right in mechanisms[left_index + 1 :]:
            output.append(
                {
                    "context": context_key,
                    "mechanism_a": left,
                    "mechanism_b": right,
                    "shared_contexts": len(contexts_by_mechanism[left] & contexts_by_mechanism[right]),
                }
            )
    output.sort(key=lambda row: (-int(row["shared_contexts"]), str(row["mechanism_a"]), str(row["mechanism_b"])))
    return output


def matched_context_effects(rows: list[dict[str, object]]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    by_context: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        by_context[str(row["cell_antigen_context"])].append(row)

    effects: list[dict[str, object]] = []
    for context, context_rows in by_context.items():
        mechanisms = sorted({str(row["payload_mechanism"]) for row in context_rows})
        if len(mechanisms) <= 1:
            continue
        by_mechanism = {
            mechanism: [row for row in context_rows if row["payload_mechanism"] == mechanism]
            for mechanism in mechanisms
        }
        for left_index, left in enumerate(mechanisms):
            for right in mechanisms[left_index + 1 :]:
                left_rows = by_mechanism[left]
                right_rows = by_mechanism[right]
                left_prev = prevalence(left_rows)
                right_prev = prevalence(right_rows)
                effects.append(
                    {
                        "cell_antigen_context": context,
                        "mechanism_a": left,
                        "mechanism_b": right,
                        "n_a": len(left_rows),
                        "n_b": len(right_rows),
                        "positives_a": sum(int(row["label_10nM"]) for row in left_rows),
                        "positives_b": sum(int(row["label_10nM"]) for row in right_rows),
                        "prevalence_a": f"{left_prev:.4f}",
                        "prevalence_b": f"{right_prev:.4f}",
                        "prevalence_diff_a_minus_b": f"{left_prev - right_prev:.4f}",
                        "payloads_a": ";".join(sorted({str(row["payload_identity"]) for row in left_rows})),
                        "payloads_b": ";".join(sorted({str(row["payload_identity"]) for row in right_rows})),
                        "has_label_variation": "1"
                        if len({str(row["label_10nM"]) for row in left_rows + right_rows}) > 1
                        else "0",
                    }
                )

    summary: list[dict[str, object]] = []
    by_pair: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in effects:
        by_pair[(str(row["mechanism_a"]), str(row["mechanism_b"]))].append(row)
    for (left, right), pair_rows in sorted(by_pair.items()):
        diffs = [float(row["prevalence_diff_a_minus_b"]) for row in pair_rows]
        weights = [int(row["n_a"]) + int(row["n_b"]) for row in pair_rows]
        weighted = sum(diff * weight for diff, weight in zip(diffs, weights)) / sum(weights)
        summary.append(
            {
                "mechanism_a": left,
                "mechanism_b": right,
                "crossed_contexts": len(pair_rows),
                "contexts_with_label_variation": sum(1 for row in pair_rows if row["has_label_variation"] == "1"),
                "observations_in_pair_contexts": sum(weights),
                "mean_prevalence_diff_a_minus_b": f"{statistics.mean(diffs):.4f}",
                "weighted_prevalence_diff_a_minus_b": f"{weighted:.4f}",
                "a_higher_contexts": sum(1 for diff in diffs if diff > 0),
                "b_higher_contexts": sum(1 for diff in diffs if diff < 0),
                "equal_contexts": sum(1 for diff in diffs if diff == 0),
            }
        )
    summary.sort(key=lambda row: (-int(row["crossed_contexts"]), str(row["mechanism_a"]), str(row["mechanism_b"])))
    effects.sort(
        key=lambda row: (
            str(row["mechanism_a"]),
            str(row["mechanism_b"]),
            -abs(float(row["prevalence_diff_a_minus_b"])),
            str(row["cell_antigen_context"]),
        )
    )
    return summary, effects


def add_context_fields(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    output = []
    for row in rows:
        enriched = dict(row)
        enriched["cell_antigen_context"] = f"{row['cell_line']} || {row['antigen_name']}"
        output.append(enriched)
    return output


@dataclass
class CategoricalNaiveBayes:
    features: list[str]
    alpha: float = 1.0

    def fit(self, rows: list[dict[str, object]]) -> None:
        self.class_counts = Counter(int(row["label_10nM"]) for row in rows)
        self.total = len(rows)
        self.values: dict[str, set[str]] = {
            feature: {str(row.get(feature, "")) for row in rows} for feature in self.features
        }
        self.feature_counts: dict[tuple[str, int, str], int] = Counter()
        for row in rows:
            y = int(row["label_10nM"])
            for feature in self.features:
                self.feature_counts[(feature, y, str(row.get(feature, "")))] += 1

    def predict_one(self, row: dict[str, object]) -> float:
        log_probs: dict[int, float] = {}
        for y in (0, 1):
            prior = (self.class_counts[y] + self.alpha) / (self.total + 2 * self.alpha)
            log_prob = math.log(prior)
            for feature in self.features:
                value = str(row.get(feature, ""))
                k = len(self.values[feature]) + 1
                count = self.feature_counts.get((feature, y, value), 0)
                prob = (count + self.alpha) / (self.class_counts[y] + self.alpha * k)
                log_prob += math.log(prob)
            log_probs[y] = log_prob
        max_log = max(log_probs.values())
        p0 = math.exp(log_probs[0] - max_log)
        p1 = math.exp(log_probs[1] - max_log)
        return p1 / (p0 + p1)

    def predict(self, rows: list[dict[str, object]]) -> list[float]:
        return [self.predict_one(row) for row in rows]


def make_random_folds(rows: list[dict[str, object]], n_folds: int, seed: int) -> list[list[int]]:
    rng = random.Random(seed)
    by_label: dict[int, list[int]] = defaultdict(list)
    for index, row in enumerate(rows):
        by_label[int(row["label_10nM"])].append(index)
    folds = [[] for _ in range(n_folds)]
    for indices in by_label.values():
        rng.shuffle(indices)
        for offset, index in enumerate(indices):
            folds[offset % n_folds].append(index)
    return folds


def make_group_folds(rows: list[dict[str, object]], group_key: str, n_folds: int, seed: int) -> list[list[int]]:
    rng = random.Random(seed)
    indices_by_group: dict[str, list[int]] = defaultdict(list)
    for index, row in enumerate(rows):
        indices_by_group[str(row[group_key])].append(index)

    groups = []
    for group, indices in indices_by_group.items():
        positives = sum(int(rows[index]["label_10nM"]) for index in indices)
        mechanisms = {str(rows[index]["payload_mechanism"]) for index in indices}
        groups.append(
            {
                "group": group,
                "indices": indices,
                "n": len(indices),
                "positives": positives,
                "mechanisms": mechanisms,
                "jitter": rng.random(),
            }
        )
    groups.sort(key=lambda item: (-int(item["n"]), -int(item["positives"]), float(item["jitter"])))

    folds = [[] for _ in range(n_folds)]
    fold_n = [0 for _ in range(n_folds)]
    fold_pos = [0 for _ in range(n_folds)]
    target_prev = sum(int(row["label_10nM"]) for row in rows) / len(rows)
    for group in groups:
        best_fold = min(
            range(n_folds),
            key=lambda fold: (
                fold_n[fold] + int(group["n"]),
                abs((fold_pos[fold] + int(group["positives"])) / max(1, fold_n[fold] + int(group["n"])) - target_prev),
                fold_n[fold],
            ),
        )
        folds[best_fold].extend(group["indices"])
        fold_n[best_fold] += int(group["n"])
        fold_pos[best_fold] += int(group["positives"])
    return folds


def rank_auc(y_true: list[int], scores: list[float]) -> float | None:
    positives = sum(y_true)
    negatives = len(y_true) - positives
    if positives == 0 or negatives == 0:
        return None

    order = sorted(range(len(scores)), key=lambda index: scores[index])
    ranks = [0.0] * len(scores)
    i = 0
    while i < len(order):
        j = i + 1
        while j < len(order) and scores[order[j]] == scores[order[i]]:
            j += 1
        avg_rank = (i + 1 + j) / 2.0
        for k in range(i, j):
            ranks[order[k]] = avg_rank
        i = j
    sum_pos_ranks = sum(ranks[index] for index, y in enumerate(y_true) if y == 1)
    return (sum_pos_ranks - positives * (positives + 1) / 2.0) / (positives * negatives)


def average_precision(y_true: list[int], scores: list[float]) -> float | None:
    positives = sum(y_true)
    if positives == 0:
        return None
    ordered = sorted(zip(scores, y_true), key=lambda item: item[0], reverse=True)
    hit_count = 0
    precision_sum = 0.0
    for rank, (_, y) in enumerate(ordered, start=1):
        if y == 1:
            hit_count += 1
            precision_sum += hit_count / rank
    return precision_sum / positives


def confusion(y_true: list[int], y_pred: list[int]) -> tuple[int, int, int, int]:
    tp = sum(1 for y, p in zip(y_true, y_pred) if y == 1 and p == 1)
    tn = sum(1 for y, p in zip(y_true, y_pred) if y == 0 and p == 0)
    fp = sum(1 for y, p in zip(y_true, y_pred) if y == 0 and p == 1)
    fn = sum(1 for y, p in zip(y_true, y_pred) if y == 1 and p == 0)
    return tp, tn, fp, fn


def balanced_accuracy(y_true: list[int], y_pred: list[int]) -> float | None:
    tp, tn, fp, fn = confusion(y_true, y_pred)
    tpr = tp / (tp + fn) if tp + fn else None
    tnr = tn / (tn + fp) if tn + fp else None
    if tpr is None or tnr is None:
        return None
    return (tpr + tnr) / 2.0


def mcc(y_true: list[int], y_pred: list[int]) -> float | None:
    tp, tn, fp, fn = confusion(y_true, y_pred)
    denom = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    if denom == 0:
        return None
    return (tp * tn - fp * fn) / denom


def log_loss(y_true: list[int], scores: list[float]) -> float:
    eps = 1e-15
    total = 0.0
    for y, score in zip(y_true, scores):
        score = min(1.0 - eps, max(eps, score))
        total -= y * math.log(score) + (1 - y) * math.log(1 - score)
    return total / len(y_true)


def brier(y_true: list[int], scores: list[float]) -> float:
    return sum((score - y) ** 2 for y, score in zip(y_true, scores)) / len(y_true)


def choose_threshold(y_true: list[int], scores: list[float]) -> float:
    candidates = sorted(set(scores))
    if not candidates:
        return 0.5
    if len(candidates) == 1:
        return candidates[0]
    best_threshold = candidates[0]
    best_key = (-999.0, -999.0)
    for threshold in candidates:
        pred = [1 if score >= threshold else 0 for score in scores]
        score_mcc = mcc(y_true, pred)
        score_ba = balanced_accuracy(y_true, pred)
        key = (
            -999.0 if score_mcc is None else score_mcc,
            -999.0 if score_ba is None else score_ba,
        )
        if key > best_key:
            best_key = key
            best_threshold = threshold
    return best_threshold


def fmt_metric(value: float | None) -> str:
    return "" if value is None or not math.isfinite(value) else f"{value:.4f}"


def evaluate_predictions(y_true: list[int], scores: list[float], threshold: float) -> dict[str, object]:
    y_pred = [1 if score >= threshold else 0 for score in scores]
    tp, tn, fp, fn = confusion(y_true, y_pred)
    return {
        "n": len(y_true),
        "positives": sum(y_true),
        "prevalence": f"{sum(y_true) / len(y_true):.4f}" if y_true else "",
        "average_precision": fmt_metric(average_precision(y_true, scores)),
        "roc_auc": fmt_metric(rank_auc(y_true, scores)),
        "balanced_accuracy": fmt_metric(balanced_accuracy(y_true, y_pred)),
        "mcc": fmt_metric(mcc(y_true, y_pred)),
        "log_loss": f"{log_loss(y_true, scores):.4f}" if y_true else "",
        "brier": f"{brier(y_true, scores):.4f}" if y_true else "",
        "threshold": f"{threshold:.6g}",
        "tp": tp,
        "tn": tn,
        "fp": fp,
        "fn": fn,
    }


def mean_metric(rows: list[dict[str, object]], key: str) -> str:
    values = []
    for row in rows:
        raw = row.get(key, "")
        if raw == "":
            continue
        try:
            values.append(float(raw))
        except (TypeError, ValueError):
            continue
    return "" if not values else f"{statistics.mean(values):.4f}"


def run_cv(
    rows: list[dict[str, object]],
    split_name: str,
    folds: list[list[int]],
    feature_sets: dict[str, list[str]],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    fold_rows: list[dict[str, object]] = []
    predictions: list[dict[str, object]] = []
    all_indices = set(range(len(rows)))
    for fold_index, test_indices in enumerate(folds):
        test_set = set(test_indices)
        train_indices = sorted(all_indices - test_set)
        train = [rows[index] for index in train_indices]
        test = [rows[index] for index in test_indices]
        y_train = [int(row["label_10nM"]) for row in train]
        y_test = [int(row["label_10nM"]) for row in test]
        for model_name, features in feature_sets.items():
            model = CategoricalNaiveBayes(features=features)
            model.fit(train)
            train_scores = model.predict(train)
            test_scores = model.predict(test)
            threshold = choose_threshold(y_train, train_scores)
            metrics = evaluate_predictions(y_test, test_scores, threshold)
            fold_rows.append(
                {
                    "split": split_name,
                    "fold": fold_index,
                    "model": model_name,
                    "features": "+".join(features) if features else "<global>",
                    **metrics,
                }
            )
            for row, score in zip(test, test_scores):
                predicted_label = 1 if score >= threshold else 0
                predictions.append(
                    {
                        "split": split_name,
                        "fold": fold_index,
                        "model": model_name,
                        "score": f"{score:.8f}",
                        "threshold": f"{threshold:.8f}",
                        "predicted_label": str(predicted_label),
                        "label_10nM": row["label_10nM"],
                        "adc_id": row["adc_id"],
                        "payload_identity": row["payload_identity"],
                        "payload_mechanism": row["payload_mechanism"],
                        "antigen_name": row["antigen_name"],
                        "cell_line": row["cell_line"],
                        "cell_antigen_context": row["cell_antigen_context"],
                    }
                )

    aggregate: list[dict[str, object]] = []
    for model_name in feature_sets:
        subset = [row for row in fold_rows if row["model"] == model_name]
        aggregate.append(
            {
                "split": split_name,
                "model": model_name,
                "features": "+".join(feature_sets[model_name]) if feature_sets[model_name] else "<global>",
                "folds": len(subset),
                "mean_average_precision": mean_metric(subset, "average_precision"),
                "mean_roc_auc": mean_metric(subset, "roc_auc"),
                "mean_balanced_accuracy": mean_metric(subset, "balanced_accuracy"),
                "mean_mcc": mean_metric(subset, "mcc"),
                "mean_log_loss": mean_metric(subset, "log_loss"),
                "mean_brier": mean_metric(subset, "brier"),
            }
        )
    return aggregate, predictions


def pooled_cv_results(predictions: list[dict[str, object]], feature_sets: dict[str, list[str]]) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    split_model_pairs = sorted({(str(row["split"]), str(row["model"])) for row in predictions})
    for split, model in split_model_pairs:
        subset = [row for row in predictions if row["split"] == split and row["model"] == model]
        y_true = [int(row["label_10nM"]) for row in subset]
        scores = [float(row["score"]) for row in subset]
        y_pred = [int(row["predicted_label"]) for row in subset]
        tp, tn, fp, fn = confusion(y_true, y_pred)
        output.append(
            {
                "split": split,
                "model": model,
                "features": "+".join(feature_sets.get(model, [])) if feature_sets.get(model, []) else "<global>",
                "aggregation": "pooled_predictions",
                "n": len(subset),
                "positives": sum(y_true),
                "prevalence": f"{sum(y_true) / len(y_true):.4f}" if y_true else "",
                "average_precision": fmt_metric(average_precision(y_true, scores)),
                "roc_auc": fmt_metric(rank_auc(y_true, scores)),
                "balanced_accuracy": fmt_metric(balanced_accuracy(y_true, y_pred)),
                "mcc": fmt_metric(mcc(y_true, y_pred)),
                "log_loss": f"{log_loss(y_true, scores):.4f}" if y_true else "",
                "brier": f"{brier(y_true, scores):.4f}" if y_true else "",
                "tp": tp,
                "tn": tn,
                "fp": fp,
                "fn": fn,
            }
        )
    return output


def fold_diagnostics(rows: list[dict[str, object]], split_name: str, folds: list[list[int]], group_key: str | None) -> list[dict[str, object]]:
    output = []
    for fold_index, indices in enumerate(folds):
        subset = [rows[index] for index in indices]
        output.append(
            {
                "split": split_name,
                "fold": fold_index,
                "n": len(subset),
                "positives": sum(int(row["label_10nM"]) for row in subset),
                "prevalence": f"{prevalence(subset):.4f}",
                "unique_payloads": unique_count(subset, "payload_identity"),
                "unique_mechanisms": unique_count(subset, "payload_mechanism"),
                "unique_antigens": unique_count(subset, "antigen_name"),
                "unique_cell_lines": unique_count(subset, "cell_line"),
                "group_key": group_key or "<row>",
                "unique_groups": "" if group_key is None else unique_count(subset, group_key),
            }
        )
    return output


def make_leave_one_group_folds(rows: list[dict[str, object]], group_key: str) -> list[list[int]]:
    indices_by_group: dict[str, list[int]] = defaultdict(list)
    for index, row in enumerate(rows):
        indices_by_group[str(row[group_key])].append(index)
    return [indices_by_group[group] for group in sorted(indices_by_group)]


def shuffled_mechanism_rows(rows: list[dict[str, object]], seed: int) -> list[dict[str, object]]:
    rng = random.Random(seed)
    mechanism_by_payload: dict[str, str] = {}
    payloads = sorted({str(row["payload_identity"]) for row in rows})
    labels = []
    for payload in payloads:
        mechanisms = {str(row["payload_mechanism"]) for row in rows if row["payload_identity"] == payload}
        labels.append(sorted(mechanisms)[0])
    shuffled = labels[:]
    rng.shuffle(shuffled)
    mechanism_by_payload = dict(zip(payloads, shuffled))
    output = []
    for row in rows:
        copied = dict(row)
        copied["payload_mechanism_shuffled"] = mechanism_by_payload[str(row["payload_identity"])]
        output.append(copied)
    return output


def run_shuffle_control(
    rows: list[dict[str, object]],
    n_folds: int,
    seed: int,
    repeats: int,
) -> list[dict[str, object]]:
    output = []
    for repeat in range(repeats):
        shuffled = shuffled_mechanism_rows(rows, seed + 1000 + repeat)
        folds = make_group_folds(shuffled, "payload_identity", n_folds, seed)
        feature_sets = {
            "bio_plus_shuffled_mechanism": [
                "antigen_name",
                "cell_line",
                "payload_mechanism_shuffled",
            ]
        }
        _, predictions = run_cv(shuffled, "payload_group", folds, feature_sets)
        row = pooled_cv_results(predictions, feature_sets)[0]
        row["repeat"] = repeat
        output.append(row)
    return output


def choose_phase2_decision(rows: list[dict[str, object]], crossed_context_rows: list[dict[str, object]]) -> tuple[str, list[str]]:
    per_mech = per_mechanism_table(rows)
    major_rows = [row for row in rows if row["payload_mechanism"] in MAJOR_MECHANISMS]
    major_per_mech = per_mechanism_table(major_rows)
    reasons = []

    min_payloads = min(int(row["unique_payloads"]) for row in major_per_mech)
    min_obs = min(int(row["observations"]) for row in major_per_mech)
    crossed_contexts = len(crossed_context_rows)
    crossed_with_label_variation = sum(1 for row in crossed_context_rows if row["has_label_variation"] == "1")
    rows_in_crossed = sum(int(row["observations"]) for row in crossed_context_rows)

    if min_payloads >= 5:
        reasons.append(f"major mechanisms each have at least 5 payloads (minimum {min_payloads}).")
    else:
        reasons.append(f"one or more major mechanisms have too few payloads (minimum {min_payloads}).")

    if min_obs >= 50:
        reasons.append(f"major mechanisms each have at least 50 observations (minimum {min_obs}).")
    else:
        reasons.append(f"one or more major mechanisms have sparse observations (minimum {min_obs}).")

    if crossed_contexts >= 20 and rows_in_crossed >= 100:
        reasons.append(
            f"there are {crossed_contexts} crossed cell-line/antigen contexts covering {rows_in_crossed} observations."
        )
    else:
        reasons.append(
            f"crossed cell-line/antigen support is limited ({crossed_contexts} contexts, {rows_in_crossed} observations)."
        )

    if crossed_with_label_variation >= 5:
        reasons.append(f"{crossed_with_label_variation} crossed contexts contain both response classes.")
    else:
        reasons.append(f"only {crossed_with_label_variation} crossed contexts contain both response classes.")

    payload_cramers = cramers_v(rows, "payload_mechanism", "payload_identity")
    if payload_cramers > 0.95:
        reasons.append(
            "payload identity almost perfectly determines mechanism; random-split mechanism gains would not be identifiable."
        )

    if min_payloads >= 5 and min_obs >= 50 and crossed_contexts >= 20 and rows_in_crossed >= 100:
        return "GO_FOR_CHEAP_FALSIFICATION", reasons
    if min_payloads >= 3 and crossed_contexts > 0:
        return "PIVOT", reasons
    return "KILL", reasons


def final_decision_from_models(model_rows: list[dict[str, object]], shuffle_rows: list[dict[str, object]]) -> tuple[str, list[str]]:
    by_split_model = {(row["split"], row["model"]): row for row in model_rows}
    reasons = []
    central_split = "payload_group"
    bio = by_split_model.get((central_split, "bio"))
    bio_mech = by_split_model.get((central_split, "bio_plus_mechanism"))
    bio_payload = by_split_model.get((central_split, "bio_plus_payload"))
    bio_payload_mech = by_split_model.get((central_split, "bio_plus_payload_plus_mechanism"))
    if not (bio and bio_mech and bio_payload and bio_payload_mech):
        return "PIVOT", ["central payload-group comparison was unavailable."]

    def metric(row: dict[str, object], key: str) -> float:
        try:
            if key in row:
                return float(row[key])
            return float(row[f"mean_{key}"])
        except (TypeError, ValueError):
            return float("nan")

    ap_gain_bio = metric(bio_mech, "average_precision") - metric(bio, "average_precision")
    ap_gain_payload = metric(bio_payload_mech, "average_precision") - metric(bio_payload, "average_precision")
    mcc_gain_bio = metric(bio_mech, "mcc") - metric(bio, "mcc")
    reasons.append(f"payload-group AP gain for bio+mechanism over bio: {ap_gain_bio:.4f}.")
    reasons.append(
        f"payload-group AP gain for bio+payload+mechanism over bio+payload: {ap_gain_payload:.4f}."
    )
    reasons.append(f"payload-group MCC gain for bio+mechanism over bio: {mcc_gain_bio:.4f}.")

    shuffled_aps = [
        metric(row, "average_precision")
        for row in shuffle_rows
        if row.get("average_precision") not in {"", None}
    ]
    if shuffled_aps:
        true_ap = metric(bio_mech, "average_precision")
        shuffled_mean = statistics.mean(shuffled_aps)
        shuffled_sorted = sorted(shuffled_aps)
        p95_index = min(len(shuffled_sorted) - 1, int(0.95 * (len(shuffled_sorted) - 1)))
        shuffled_p95 = shuffled_sorted[p95_index]
        reasons.append(
            f"true bio+mechanism AP {true_ap:.4f}; shuffled-mechanism mean {shuffled_mean:.4f}, p95 {shuffled_p95:.4f}."
        )
        if true_ap <= shuffled_p95:
            reasons.append("true mechanism does not beat the shuffled-mechanism control.")
            return "PIVOT", reasons

    if ap_gain_bio >= 0.02 and ap_gain_payload >= 0.02 and mcc_gain_bio > 0:
        return "GO", reasons
    return "PIVOT", reasons


def svg_bar_chart(path: Path, rows: list[dict[str, object]], label_key: str, value_key: str, title: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    width = 820
    bar_h = 32
    gap = 12
    left = 260
    right = 40
    top = 56
    height = top + len(rows) * (bar_h + gap) + 28
    values = [float(row[value_key]) for row in rows]
    max_value = max(values) if values else 1.0
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        f'<text x="24" y="32" font-family="Arial, sans-serif" font-size="20" font-weight="700" fill="#1f2937">{title}</text>',
    ]
    palette = ["#2f6f73", "#b84a62", "#7a5c1f", "#4f67a5", "#6b7280"]
    for idx, row in enumerate(rows):
        y = top + idx * (bar_h + gap)
        value = float(row[value_key])
        bar_w = (width - left - right) * (value / max_value if max_value else 0)
        color = palette[idx % len(palette)]
        label = str(row[label_key]).replace("&", "&amp;")
        lines.append(f'<text x="24" y="{y + 22}" font-family="Arial, sans-serif" font-size="14" fill="#111827">{label}</text>')
        lines.append(f'<rect x="{left}" y="{y}" width="{bar_w:.1f}" height="{bar_h}" fill="{color}"/>')
        lines.append(f'<text x="{left + bar_w + 8}" y="{y + 22}" font-family="Arial, sans-serif" font-size="14" fill="#111827">{value:.3g}</text>')
    lines.append("</svg>")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_report(
    path: Path,
    raw_summary: list[dict[str, object]],
    xlsx_summary: list[dict[str, object]],
    activity_rows: list[dict[str, str]],
    usable_rows: list[dict[str, object]],
    exclusions: list[dict[str, object]],
    collapsed: list[dict[str, object]],
    conflicts: list[dict[str, object]],
    per_mech: list[dict[str, object]],
    association: list[dict[str, object]],
    overlap_summaries: dict[str, list[dict[str, object]]],
    overlap_examples: dict[str, list[dict[str, object]]],
    pairwise_overlaps: list[dict[str, object]],
    matched_summary: list[dict[str, object]],
    matched_effects: list[dict[str, object]],
    phase2_decision: str,
    phase2_reasons: list[str],
    model_results: list[dict[str, object]],
    shuffle_results: list[dict[str, object]],
    final_decision: str,
    final_reasons: list[str],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    response_counts = Counter(str(row["label_10nM"]) for row in collapsed)
    exclusion_counts = Counter(str(row["reason"]) for row in exclusions)

    lines = [
        "# ADC Mechanism Identifiability Study",
        "",
        "## Scope",
        "",
        "- Raw data were read from local `data/raw/`; no replacement dataset was downloaded.",
        "- No neural model or complex architecture was trained.",
        f"- Parsed ADCdb cell-line activity rows available: {len(activity_rows):,}.",
        f"- Usable deterministic molar potency rows at 10 nM: {len(usable_rows):,}.",
        f"- Collapsed modeling observations: {len(collapsed):,}.",
        f"- Collapsed response distribution: positives={response_counts['1']:,}, negatives={response_counts['0']:,}, prevalence={prevalence(collapsed):.3f}.",
        f"- Duplicate collapse conflicts at 10 nM: {len(conflicts):,} context groups.",
        "",
        "Usable rows require potency endpoint (IC50/LC50/EC50/GI50), molar unit convertible to nM, and a determinate 10 nM label. Mass units were excluded because molecular-weight conversion is not available in these local rows.",
        "",
        "## Raw Data Inventory",
        "",
        markdown_table(raw_summary, ["top_level", "files", "total_bytes", "suffixes"]),
        "",
        "## ADCNet Excel Files",
        "",
        markdown_table(xlsx_summary, ["path", "rows", "columns", "label_columns", "label_counts"]),
        "",
        "Local interpretation: `data.xlsx` and `files_data.xlsx` are byte-identical ADC-level ADCNet tables. They include threshold labels but no cell-line activity rows. `t_data.xlsx` is smaller and lacks a 10 nM label, so it is not the primary dataset for this question.",
        "",
        "## Endpoint And Label Accounting",
        "",
        markdown_table(
            [{"reason": key, "rows": exclusion_counts[key]} for key in sorted(exclusion_counts)],
            ["reason",
             "rows"],
        ),
        "",
        "## Per-Mechanism Audit",
        "",
        markdown_table(
            per_mech,
            [
                "payload_mechanism",
                "observations",
                "positives_10nM",
                "prevalence_10nM",
                "unique_payloads",
                "effective_payload_n",
                "unique_antigens",
                "unique_cell_lines",
                "top_payload",
                "top_payload_share",
                "top_antigen",
                "top_antigen_share",
            ],
        ),
        "",
        "## Mechanism Associations",
        "",
        markdown_table(association, ["association", "right_unique_values", "cramers_v"]),
        "",
        "Cramer's V near 1.0 indicates strong association. Payload identity almost deterministically implies mechanism, so seen-payload random splits cannot establish mechanism-level generalization.",
        "",
        "## Cross-Mechanism Biological Overlap",
        "",
    ]
    for key, summary in overlap_summaries.items():
        lines.extend(
            [
                f"### {key}",
                "",
                markdown_table(summary, ["context", "mechanism_count", "unique_contexts"]),
                "",
            ]
        )
    lines.extend(
        [
            "Pairwise overlaps:",
            "",
            markdown_table(pairwise_overlaps, ["context", "mechanism_a", "mechanism_b", "shared_contexts"]),
            "",
            "Top crossed cell-line/antigen contexts:",
            "",
            markdown_table(
                overlap_examples["cell_antigen_context"][:20],
                [
                    "context",
                    "mechanism_count",
                    "observations",
                    "unique_payloads",
                    "has_label_variation",
                    "mechanisms",
                ],
            ),
            "",
            "## Matched-Context Mechanism Effects",
            "",
            "This uses only observations where the same cell-line/antigen context has more than one major mechanism. It is a low-power stress test for mechanism effects after holding both measured cell line and ADC antigen fixed.",
            "",
            markdown_table(
                matched_summary,
                [
                    "mechanism_a",
                    "mechanism_b",
                    "crossed_contexts",
                    "contexts_with_label_variation",
                    "observations_in_pair_contexts",
                    "mean_prevalence_diff_a_minus_b",
                    "weighted_prevalence_diff_a_minus_b",
                    "a_higher_contexts",
                    "b_higher_contexts",
                    "equal_contexts",
                ],
            ),
            "",
            "Largest matched-context differences:",
            "",
            markdown_table(
                matched_effects[:15],
                [
                    "cell_antigen_context",
                    "mechanism_a",
                    "mechanism_b",
                    "n_a",
                    "n_b",
                    "prevalence_a",
                    "prevalence_b",
                    "prevalence_diff_a_minus_b",
                    "payloads_a",
                    "payloads_b",
                ],
            ),
            "",
            "## Phase 2 Decision",
            "",
            f"- Decision before modeling: **{phase2_decision}**",
            "",
        ]
    )
    for reason in phase2_reasons:
        lines.append(f"- {reason}")
    lines.extend(
        [
            "",
            "## Cheap Modeling Falsification",
            "",
            "Model: categorical Naive Bayes with Laplace smoothing. Splits are deterministic. The central split is grouped by payload identity, so test payloads are unseen during training.",
            "",
            markdown_table(
                model_results,
                [
                    "split",
                    "model",
                    "average_precision",
                    "roc_auc",
                    "balanced_accuracy",
                    "mcc",
                    "log_loss",
                    "brier",
                ],
            ),
            "",
            "Shuffle control: mechanism labels were permuted at payload level and evaluated under the same payload-group split.",
            "",
        ]
    )
    if shuffle_results:
        shuffled_aps = [float(row["average_precision"]) for row in shuffle_results if row["average_precision"]]
        lines.extend(
            [
                f"- Shuffle repeats: {len(shuffle_results)}",
                f"- Shuffled AP mean: {statistics.mean(shuffled_aps):.4f}",
                f"- Shuffled AP min/max: {min(shuffled_aps):.4f} / {max(shuffled_aps):.4f}",
                "",
            ]
        )
    lines.extend(
        [
            "## Current Decision",
            "",
            f"- **{final_decision}**",
            "",
        ]
    )
    for reason in final_reasons:
        lines.append(f"- {reason}")
    lines.extend(
        [
            "",
            "## Pivot Recommendation",
            "",
            "The defensible reformulation is narrower: test whether coarse payload mechanism can act as a low-dimensional prior for unseen payloads in ADCdb-like metadata, not whether it independently explains ADC response beyond payload identity and biological context.",
            "",
            "## Highest-Value Next Step",
            "",
            "Stop model escalation. The next useful step is curation, not architecture: validate cell-line naming/synonyms, endpoint aggregation, and mass-unit conversion assumptions to see whether more genuinely matched cell-line/antigen contexts can be recovered from the existing local ADCdb pages.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    raw_inventory, raw_summary = summarize_raw_files(args.raw_dir)
    xlsx_summary = summarize_adcnet_xlsx(args.raw_dir)
    write_csv(args.out_dir / "raw_data_inventory.csv", raw_inventory)
    write_csv(args.out_dir / "raw_data_summary.csv", raw_summary)
    write_csv(args.out_dir / "adcnet_excel_summary.csv", xlsx_summary)

    activity_rows = read_csv(args.activity_csv)
    usable_rows, exclusions = derive_usable_rows(activity_rows)
    collapsed, conflicts = collapse_context_rows(usable_rows)
    collapsed = add_context_fields(collapsed)
    major_rows = [row for row in collapsed if row["payload_mechanism"] in MAJOR_MECHANISMS]

    context_fieldnames = [
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
        "cell_line",
        "disease_model",
        "value_nM_min",
        "label_10nM",
        "source_activity_rows",
        "has_conflicting_10nM_labels",
        "source_activity_ids",
        "source_url",
        "cell_antigen_context",
    ]
    write_csv(args.context_csv, collapsed, context_fieldnames)
    write_csv(args.out_dir / "modeling_dataset_10nM_context.csv", collapsed, context_fieldnames)
    write_csv(args.out_dir / "excluded_activity_rows_10nM.csv", exclusions)
    write_csv(args.out_dir / "duplicate_label_conflicts_10nM.csv", conflicts)

    per_mech = per_mechanism_table(collapsed)
    major_per_mech = per_mechanism_table(major_rows)
    association = association_table(collapsed)
    write_csv(args.out_dir / "mechanism_summary_10nM.csv", per_mech)
    write_csv(args.out_dir / "major_mechanism_summary_10nM.csv", major_per_mech)
    write_csv(args.out_dir / "mechanism_associations.csv", association)
    write_csv(args.out_dir / "mechanism_payload_matrix.csv", matrix_counts(collapsed, "payload_mechanism", "payload_identity"))
    write_csv(args.out_dir / "mechanism_antigen_matrix.csv", matrix_counts(collapsed, "payload_mechanism", "antigen_name"))

    overlap_summaries: dict[str, list[dict[str, object]]] = {}
    overlap_examples: dict[str, list[dict[str, object]]] = {}
    pairwise: list[dict[str, object]] = []
    for context_key in ["cell_line", "antigen_name", "cell_antigen_context"]:
        summary, examples = overlap_table(collapsed, context_key)
        overlap_summaries[context_key] = summary
        overlap_examples[context_key] = examples
        pairwise.extend(pairwise_context_overlap(collapsed, context_key))
        write_csv(args.out_dir / f"{context_key}_mechanism_overlap_summary.csv", summary)
        write_csv(args.out_dir / f"{context_key}_crossed_examples.csv", examples)
    write_csv(args.out_dir / "pairwise_context_overlaps.csv", pairwise)

    _, major_crossed_contexts = overlap_table(major_rows, "cell_antigen_context")
    phase2_decision, phase2_reasons = choose_phase2_decision(major_rows, major_crossed_contexts)
    matched_summary, matched_effects = matched_context_effects(major_rows)
    write_csv(args.out_dir / "matched_context_mechanism_effect_summary.csv", matched_summary)
    write_csv(args.out_dir / "matched_context_mechanism_effects.csv", matched_effects)

    feature_sets = {
        "global": [],
        "mechanism_only": ["payload_mechanism"],
        "payload_only": ["payload_identity"],
        "antigen_only": ["antigen_name"],
        "cell_only": ["cell_line"],
        "bio": ["antigen_name", "cell_line"],
        "bio_plus_mechanism": ["antigen_name", "cell_line", "payload_mechanism"],
        "bio_plus_payload": ["antigen_name", "cell_line", "payload_identity"],
        "bio_plus_payload_plus_mechanism": [
            "antigen_name",
            "cell_line",
            "payload_identity",
            "payload_mechanism",
        ],
    }

    model_rows: list[dict[str, object]] = []
    prediction_rows: list[dict[str, object]] = []
    fold_rows: list[dict[str, object]] = []
    if phase2_decision == "GO_FOR_CHEAP_FALSIFICATION":
        split_defs = {
            "random_row": (make_random_folds(major_rows, args.folds, args.seed), None),
            "payload_group": (
                make_group_folds(major_rows, "payload_identity", args.folds, args.seed),
                "payload_identity",
            ),
            "payload_leave_one_out": (
                make_leave_one_group_folds(major_rows, "payload_identity"),
                "payload_identity",
            ),
            "cell_line_group": (make_group_folds(major_rows, "cell_line", args.folds, args.seed), "cell_line"),
            "cell_antigen_context_group": (
                make_group_folds(major_rows, "cell_antigen_context", args.folds, args.seed),
                "cell_antigen_context",
            ),
        }
        for split_name, (folds, group_key) in split_defs.items():
            aggregate, predictions = run_cv(major_rows, split_name, folds, feature_sets)
            model_rows.extend(aggregate)
            prediction_rows.extend(predictions)
            fold_rows.extend(fold_diagnostics(major_rows, split_name, folds, group_key))
        shuffle_rows = run_shuffle_control(major_rows, args.folds, args.seed, args.shuffle_repeats)
    else:
        shuffle_rows = []

    pooled_model_rows = pooled_cv_results(prediction_rows, feature_sets) if prediction_rows else []
    write_csv(args.out_dir / "cheap_model_results.csv", pooled_model_rows)
    write_csv(args.out_dir / "cheap_model_fold_mean_results.csv", model_rows)
    write_csv(args.out_dir / "cheap_model_predictions.csv", prediction_rows)
    write_csv(args.out_dir / "cv_fold_diagnostics.csv", fold_rows)
    write_csv(args.out_dir / "shuffled_mechanism_control.csv", shuffle_rows)

    central_errors = []
    for row in prediction_rows:
        if row["split"] != "payload_group" or row["model"] != "bio_plus_mechanism":
            continue
        score = float(row["score"])
        label = int(row["label_10nM"])
        predicted = int(row["predicted_label"])
        if predicted != label:
            copied = dict(row)
            copied["error_type"] = "false_positive" if predicted == 1 else "false_negative"
            copied["confidence"] = f"{abs(score - float(row['threshold'])):.8f}"
            central_errors.append(copied)
    central_errors.sort(key=lambda row: -float(row["confidence"]))
    write_csv(args.out_dir / "representative_errors_payload_group.csv", central_errors[:50])

    final_decision, final_reasons = (
        final_decision_from_models(pooled_model_rows, shuffle_rows)
        if phase2_decision == "GO_FOR_CHEAP_FALSIFICATION"
        else (phase2_decision, phase2_reasons)
    )
    if matched_summary:
        variable_contexts = sum(int(row["contexts_with_label_variation"]) for row in matched_summary)
        total_pair_contexts = sum(int(row["crossed_contexts"]) for row in matched_summary)
        final_reasons.append(
            f"matched-context stress test has {variable_contexts}/{total_pair_contexts} pairwise contexts with label variation."
        )

    svg_bar_chart(
        args.out_dir / "figures" / "mechanism_observation_counts.svg",
        per_mech,
        "payload_mechanism",
        "observations",
        "Collapsed 10 nM Observations By Mechanism",
    )
    prevalence_rows = [
        dict(row, prevalence_float=float(row["prevalence_10nM"]))
        for row in per_mech
    ]
    svg_bar_chart(
        args.out_dir / "figures" / "mechanism_response_prevalence.svg",
        prevalence_rows,
        "payload_mechanism",
        "prevalence_float",
        "10 nM Response Prevalence By Mechanism",
    )

    write_report(
        args.out_dir / "identifiability_study.md",
        raw_summary,
        xlsx_summary,
        activity_rows,
        usable_rows,
        exclusions,
        collapsed,
        conflicts,
        per_mech,
        association,
        overlap_summaries,
        overlap_examples,
        pairwise,
        matched_summary,
        matched_effects,
        phase2_decision,
        phase2_reasons,
        pooled_model_rows,
        shuffle_rows,
        final_decision,
        final_reasons,
    )

    print(f"Wrote {args.out_dir / 'identifiability_study.md'}")
    print(f"Wrote {args.context_csv}")
    print(f"Phase 2 decision: {phase2_decision}")
    print(f"Current decision: {final_decision}")


if __name__ == "__main__":
    main()
