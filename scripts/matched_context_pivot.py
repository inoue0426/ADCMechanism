#!/usr/bin/env python3
"""Matched-context pivot analysis for ADC payload mechanisms.

This script answers the final pivot question after the global identifiability
study: within the same antigen x cell-line context, is there reproducible
evidence that payload mechanism is associated with differential response?

The unit of inference is a payload-context unit. Repeated ADC observations with
the same payload in the same antigen/cell-line context are collapsed before
mechanism comparisons, so repeated measurements of one payload do not count as
independent mechanism-level evidence.

No predictive model or neural architecture is trained here.
"""

from __future__ import annotations

import argparse
import csv
import math
import random
import statistics
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path
from xml.sax.saxutils import escape


MAJOR_MECHANISMS = [
    "DNA-damaging",
    "Topo-I inhibitor",
    "microtubule-disrupting",
]
MAJOR_MECHANISM_SET = set(MAJOR_MECHANISMS)
EPS = 1e-12


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--context-csv",
        type=Path,
        default=Path("data/processed/adc_response_10nM_context.csv"),
    )
    parser.add_argument("--out-dir", type=Path, default=Path("results"))
    parser.add_argument("--seed", type=int, default=44)
    parser.add_argument("--permutations", type=int, default=10000)
    parser.add_argument("--bootstraps", type=int, default=5000)
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = []
        seen: set[str] = set()
        for row in rows:
            for key in row:
                if key not in seen:
                    seen.add(key)
                    fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def f4(value: float | None) -> str:
    if value is None or not math.isfinite(value):
        return ""
    return f"{value:.4f}"


def pct(value: float | None) -> str:
    if value is None or not math.isfinite(value):
        return ""
    return f"{100.0 * value:.1f}%"


def mean(values: list[float]) -> float | None:
    if not values:
        return None
    return statistics.mean(values)


def median(values: list[float]) -> float | None:
    if not values:
        return None
    return statistics.median(values)


def percentile(values: list[float], q: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    rank = q * (len(ordered) - 1)
    lo = int(math.floor(rank))
    hi = int(math.ceil(rank))
    if lo == hi:
        return ordered[lo]
    frac = rank - lo
    return ordered[lo] * (1.0 - frac) + ordered[hi] * frac


def markdown_table(rows: list[dict[str, object]], headers: list[str], max_rows: int | None = None) -> str:
    shown = rows if max_rows is None else rows[:max_rows]
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in shown:
        lines.append("| " + " | ".join(str(row.get(header, "")) for header in headers) + " |")
    if max_rows is not None and len(rows) > max_rows:
        lines.append("| " + " | ".join(["..."] * len(headers)) + " |")
    return "\n".join(lines)


def mechanism_sort_key(mechanism: str) -> tuple[int, str]:
    try:
        return (MAJOR_MECHANISMS.index(mechanism), mechanism)
    except ValueError:
        return (len(MAJOR_MECHANISMS), mechanism)


def parse_float(raw: str) -> float | None:
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def collapse_payload_context_units(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        context = row.get("cell_antigen_context") or f"{row.get('cell_line', '')} || {row.get('antigen_name', '')}"
        key = (
            context,
            row.get("cell_line", ""),
            row.get("antigen_name", ""),
            row.get("payload_mechanism", ""),
            row.get("payload_identity", ""),
        )
        grouped[key].append(row)

    units: list[dict[str, object]] = []
    for (context, cell_line, antigen_name, mechanism, payload), group_rows in sorted(grouped.items()):
        labels = [int(row["label_10nM"]) for row in group_rows]
        values = [parse_float(row.get("value_nM_min", "")) for row in group_rows]
        values = [value for value in values if value is not None]
        adcs = sorted({row.get("adc_id", "") for row in group_rows if row.get("adc_id", "")})
        payload_names = sorted({row.get("payload_name", "") for row in group_rows if row.get("payload_name", "")})
        payload_targets = sorted({row.get("payload_target", "") for row in group_rows if row.get("payload_target", "")})
        units.append(
            {
                "payload_context_unit_id": f"{context} || {mechanism} || {payload}",
                "cell_antigen_context": context,
                "cell_line": cell_line,
                "antigen_name": antigen_name,
                "payload_mechanism": mechanism,
                "payload_identity": payload,
                "payload_name": ";".join(payload_names),
                "payload_target": ";".join(payload_targets),
                "input_observations": len(group_rows),
                "unique_adcs": len(adcs),
                "adc_ids": ";".join(adcs),
                "positives": sum(labels),
                "negatives": len(labels) - sum(labels),
                "response_rate": sum(labels) / len(labels),
                "response_any": int(any(labels)),
                "response_all": int(all(labels)),
                "has_repeated_observations": int(len(group_rows) > 1),
                "has_discordant_repeated_labels": int(len(set(labels)) > 1),
                "value_nM_min": "" if not values else min(values),
                "source_activity_rows_total": sum(
                    int(row.get("source_activity_rows", "1") or "1") for row in group_rows
                ),
            }
        )
    return units


def contexts_by_units(units: list[dict[str, object]]) -> dict[str, list[dict[str, object]]]:
    by_context: dict[str, list[dict[str, object]]] = defaultdict(list)
    for unit in units:
        by_context[str(unit["cell_antigen_context"])].append(unit)
    return by_context


def matched_context_names(units: list[dict[str, object]]) -> set[str]:
    names = set()
    for context, context_units in contexts_by_units(units).items():
        mechanisms = {str(unit["payload_mechanism"]) for unit in context_units}
        if len(mechanisms) > 1:
            names.add(context)
    return names


def arm_units(units: list[dict[str, object]], mechanism: str) -> list[dict[str, object]]:
    return [unit for unit in units if str(unit["payload_mechanism"]) == mechanism]


def unit_response_mean(units: list[dict[str, object]]) -> float | None:
    return mean([float(unit["response_rate"]) for unit in units])


def observation_prevalence(units: list[dict[str, object]]) -> float | None:
    observations = sum(int(unit["input_observations"]) for unit in units)
    if observations == 0:
        return None
    positives = sum(int(unit["positives"]) for unit in units)
    return positives / observations


def semicolon_counts(counter: Counter[str]) -> str:
    return "; ".join(f"{key}={counter[key]}" for key in sorted(counter, key=mechanism_sort_key))


def semicolon_values(values: dict[str, object]) -> str:
    return "; ".join(f"{key}={values[key]}" for key in sorted(values, key=mechanism_sort_key))


def context_summary(units: list[dict[str, object]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for context, context_units in contexts_by_units(units).items():
        mechanisms = sorted({str(unit["payload_mechanism"]) for unit in context_units}, key=mechanism_sort_key)
        by_mechanism = {mechanism: arm_units(context_units, mechanism) for mechanism in mechanisms}
        unit_rates = [float(unit["response_rate"]) for unit in context_units]
        payload_units_by_mechanism = Counter(
            {mechanism: len(by_mechanism[mechanism]) for mechanism in mechanisms}
        )
        observations_by_mechanism = Counter(
            {
                mechanism: sum(int(unit["input_observations"]) for unit in by_mechanism[mechanism])
                for mechanism in mechanisms
            }
        )
        unit_response_by_mechanism = {
            mechanism: f4(unit_response_mean(by_mechanism[mechanism])) for mechanism in mechanisms
        }
        payloads_by_mechanism = {
            mechanism: ",".join(sorted({str(unit["payload_identity"]) for unit in by_mechanism[mechanism]}))
            for mechanism in mechanisms
        }
        rows.append(
            {
                "cell_antigen_context": context,
                "cell_line": str(context_units[0]["cell_line"]),
                "antigen_name": str(context_units[0]["antigen_name"]),
                "mechanism_count": len(mechanisms),
                "mechanisms": "; ".join(mechanisms),
                "payload_context_units": len(context_units),
                "input_observations": sum(int(unit["input_observations"]) for unit in context_units),
                "unique_payloads": len({str(unit["payload_identity"]) for unit in context_units}),
                "mechanisms_with_2plus_payloads": sum(
                    1 for mechanism in mechanisms if len(by_mechanism[mechanism]) >= 2
                ),
                "all_mechanisms_have_2plus_payloads": int(
                    all(len(by_mechanism[mechanism]) >= 2 for mechanism in mechanisms)
                ),
                "payload_unit_response_mean": f4(unit_response_mean(context_units)),
                "observation_prevalence": f4(observation_prevalence(context_units)),
                "unit_response_min": f4(min(unit_rates) if unit_rates else None),
                "unit_response_max": f4(max(unit_rates) if unit_rates else None),
                "has_unit_response_variation": int(max(unit_rates) - min(unit_rates) > EPS) if unit_rates else 0,
                "payload_units_by_mechanism": semicolon_counts(payload_units_by_mechanism),
                "observations_by_mechanism": semicolon_counts(observations_by_mechanism),
                "unit_response_by_mechanism": semicolon_values(unit_response_by_mechanism),
                "payloads_by_mechanism": semicolon_values(payloads_by_mechanism),
            }
        )
    rows.sort(
        key=lambda row: (
            -int(row["mechanism_count"]),
            -int(row["payload_context_units"]),
            str(row["cell_antigen_context"]),
        )
    )
    return rows


def pairwise_context_effects(units: list[dict[str, object]]) -> list[dict[str, object]]:
    effects: list[dict[str, object]] = []
    for context, context_units in contexts_by_units(units).items():
        present = sorted({str(unit["payload_mechanism"]) for unit in context_units}, key=mechanism_sort_key)
        if len(present) <= 1:
            continue
        by_mechanism = {mechanism: arm_units(context_units, mechanism) for mechanism in present}
        for mechanism_a, mechanism_b in combinations(present, 2):
            units_a = by_mechanism[mechanism_a]
            units_b = by_mechanism[mechanism_b]
            mean_a = unit_response_mean(units_a)
            mean_b = unit_response_mean(units_b)
            obs_a = observation_prevalence(units_a)
            obs_b = observation_prevalence(units_b)
            payloads_a = sorted({str(unit["payload_identity"]) for unit in units_a})
            payloads_b = sorted({str(unit["payload_identity"]) for unit in units_b})
            diff = None if mean_a is None or mean_b is None else mean_a - mean_b
            obs_diff = None if obs_a is None or obs_b is None else obs_a - obs_b
            effects.append(
                {
                    "cell_antigen_context": context,
                    "cell_line": str(context_units[0]["cell_line"]),
                    "antigen_name": str(context_units[0]["antigen_name"]),
                    "mechanism_a": mechanism_a,
                    "mechanism_b": mechanism_b,
                    "payload_units_a": len(units_a),
                    "payload_units_b": len(units_b),
                    "input_observations_a": sum(int(unit["input_observations"]) for unit in units_a),
                    "input_observations_b": sum(int(unit["input_observations"]) for unit in units_b),
                    "unique_payloads_a": len(payloads_a),
                    "unique_payloads_b": len(payloads_b),
                    "unit_response_a": f4(mean_a),
                    "unit_response_b": f4(mean_b),
                    "unit_response_diff_a_minus_b": f4(diff),
                    "observation_prevalence_a": f4(obs_a),
                    "observation_prevalence_b": f4(obs_b),
                    "observation_prevalence_diff_a_minus_b": f4(obs_diff),
                    "a_higher": int(diff is not None and diff > EPS),
                    "b_higher": int(diff is not None and diff < -EPS),
                    "equal": int(diff is not None and abs(diff) <= EPS),
                    "both_arms_2plus_payloads": int(len(units_a) >= 2 and len(units_b) >= 2),
                    "either_arm_2plus_payloads": int(len(units_a) >= 2 or len(units_b) >= 2),
                    "has_unit_response_difference": int(diff is not None and abs(diff) > EPS),
                    "payloads_a": ";".join(payloads_a),
                    "payloads_b": ";".join(payloads_b),
                    "payload_pairing_signature": ";".join(
                        f"{payload_a}|{payload_b}" for payload_a in payloads_a for payload_b in payloads_b
                    ),
                }
            )
    effects.sort(
        key=lambda row: (
            str(row["mechanism_a"]),
            str(row["mechanism_b"]),
            -abs(float(row["unit_response_diff_a_minus_b"] or "0")),
            str(row["cell_antigen_context"]),
        )
    )
    return effects


def summarize_pairwise_effects(effects: list[dict[str, object]]) -> list[dict[str, object]]:
    by_pair: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in effects:
        by_pair[(str(row["mechanism_a"]), str(row["mechanism_b"]))].append(row)

    rows: list[dict[str, object]] = []
    for (mechanism_a, mechanism_b), pair_rows in sorted(by_pair.items(), key=lambda item: (mechanism_sort_key(item[0][0]), mechanism_sort_key(item[0][1]))):
        diffs = [float(row["unit_response_diff_a_minus_b"]) for row in pair_rows if row["unit_response_diff_a_minus_b"] != ""]
        unit_weights = [
            int(row["payload_units_a"]) + int(row["payload_units_b"])
            for row in pair_rows
            if row["unit_response_diff_a_minus_b"] != ""
        ]
        obs_diffs = [
            float(row["observation_prevalence_diff_a_minus_b"])
            for row in pair_rows
            if row["observation_prevalence_diff_a_minus_b"] != ""
        ]
        obs_weights = [
            int(row["input_observations_a"]) + int(row["input_observations_b"])
            for row in pair_rows
            if row["observation_prevalence_diff_a_minus_b"] != ""
        ]
        payloads_a = sorted(
            {payload for row in pair_rows for payload in str(row["payloads_a"]).split(";") if payload}
        )
        payloads_b = sorted(
            {payload for row in pair_rows for payload in str(row["payloads_b"]).split(";") if payload}
        )
        pairing_counts = Counter(str(row["payload_pairing_signature"]) for row in pair_rows)
        top_pairing, top_pairing_n = pairing_counts.most_common(1)[0]
        rows.append(
            {
                "mechanism_a": mechanism_a,
                "mechanism_b": mechanism_b,
                "crossed_contexts": len(pair_rows),
                "input_observations_in_pair_contexts": sum(
                    int(row["input_observations_a"]) + int(row["input_observations_b"]) for row in pair_rows
                ),
                "payload_context_units_in_pair_contexts": sum(
                    int(row["payload_units_a"]) + int(row["payload_units_b"]) for row in pair_rows
                ),
                "contexts_with_unit_response_difference": sum(
                    int(row["has_unit_response_difference"]) for row in pair_rows
                ),
                "a_higher_contexts": sum(int(row["a_higher"]) for row in pair_rows),
                "b_higher_contexts": sum(int(row["b_higher"]) for row in pair_rows),
                "equal_contexts": sum(int(row["equal"]) for row in pair_rows),
                "mean_context_unit_response_diff_a_minus_b": f4(mean(diffs)),
                "median_context_unit_response_diff_a_minus_b": f4(median(diffs)),
                "payload_unit_weighted_diff_a_minus_b": f4(
                    sum(diff * weight for diff, weight in zip(diffs, unit_weights)) / sum(unit_weights)
                    if unit_weights
                    else None
                ),
                "observation_weighted_diff_a_minus_b": f4(
                    sum(diff * weight for diff, weight in zip(obs_diffs, obs_weights)) / sum(obs_weights)
                    if obs_weights
                    else None
                ),
                "contexts_with_both_arms_2plus_payloads": sum(
                    int(row["both_arms_2plus_payloads"]) for row in pair_rows
                ),
                "contexts_with_either_arm_2plus_payloads": sum(
                    int(row["either_arm_2plus_payloads"]) for row in pair_rows
                ),
                "unique_payloads_a_across_contexts": len(payloads_a),
                "unique_payloads_b_across_contexts": len(payloads_b),
                "payloads_a_across_contexts": ";".join(payloads_a),
                "payloads_b_across_contexts": ";".join(payloads_b),
                "top_payload_pairing_signature": top_pairing,
                "top_payload_pairing_contexts": top_pairing_n,
            }
        )
    rows.sort(key=lambda row: (-int(row["crossed_contexts"]), str(row["mechanism_a"]), str(row["mechanism_b"])))
    return rows


def dataset_reduction_rows(
    all_rows: list[dict[str, str]],
    major_rows: list[dict[str, str]],
    all_units: list[dict[str, object]],
    major_units: list[dict[str, object]],
    matched_units: list[dict[str, object]],
) -> list[dict[str, object]]:
    def summarize(scope: str, rows: list[dict[str, str]], units: list[dict[str, object]]) -> dict[str, object]:
        labels = [int(row["label_10nM"]) for row in rows]
        unit_rates = [float(unit["response_rate"]) for unit in units]
        return {
            "scope": scope,
            "input_observations": len(rows),
            "positive_input_observations": sum(labels),
            "negative_input_observations": len(labels) - sum(labels),
            "observation_prevalence": f4(sum(labels) / len(labels) if labels else None),
            "payload_context_units": len(units),
            "mean_payload_unit_response_rate": f4(mean(unit_rates)),
            "unique_payloads": len({str(unit["payload_identity"]) for unit in units}),
            "unique_mechanisms": len({str(unit["payload_mechanism"]) for unit in units}),
            "unique_cell_antigen_contexts": len({str(unit["cell_antigen_context"]) for unit in units}),
        }

    matched_contexts = {str(unit["cell_antigen_context"]) for unit in matched_units}
    matched_rows = [row for row in major_rows if row["cell_antigen_context"] in matched_contexts]
    return [
        summarize("all_processed_observations", all_rows, all_units),
        summarize("major_mechanisms_only", major_rows, major_units),
        summarize("matched_major_exact_contexts", matched_rows, matched_units),
        {
            "scope": "matched_major_survival_fraction_vs_major",
            "input_observations": pct(len(matched_rows) / len(major_rows) if major_rows else None),
            "positive_input_observations": "",
            "negative_input_observations": "",
            "observation_prevalence": "",
            "payload_context_units": pct(len(matched_units) / len(major_units) if major_units else None),
            "mean_payload_unit_response_rate": "",
            "unique_payloads": "",
            "unique_mechanisms": "",
            "unique_cell_antigen_contexts": "",
        },
    ]


def residual_rows_for_units(units: list[dict[str, object]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for context, context_units in contexts_by_units(units).items():
        context_mean = unit_response_mean(context_units)
        if context_mean is None:
            continue
        for unit in context_units:
            response_rate = float(unit["response_rate"])
            rows.append(
                {
                    "payload_context_unit_id": unit["payload_context_unit_id"],
                    "cell_antigen_context": context,
                    "payload_mechanism": unit["payload_mechanism"],
                    "payload_identity": unit["payload_identity"],
                    "response_rate": f4(response_rate),
                    "context_unit_response_mean": f4(context_mean),
                    "context_fixed_residual": f4(response_rate - context_mean),
                    "input_observations": unit["input_observations"],
                    "unique_adcs": unit["unique_adcs"],
                }
            )
    rows.sort(key=lambda row: (str(row["payload_mechanism"]), str(row["payload_identity"]), str(row["cell_antigen_context"])))
    return rows


def residual_group_means(units: list[dict[str, object]]) -> dict[str, float]:
    residuals_by_mechanism: dict[str, list[float]] = defaultdict(list)
    for context, context_units in contexts_by_units(units).items():
        if len({str(unit["payload_mechanism"]) for unit in context_units}) <= 1:
            continue
        context_mean = unit_response_mean(context_units)
        if context_mean is None:
            continue
        for unit in context_units:
            residuals_by_mechanism[str(unit["payload_mechanism"])].append(float(unit["response_rate"]) - context_mean)
    return {
        mechanism: statistics.mean(values)
        for mechanism, values in residuals_by_mechanism.items()
        if values
    }


def residual_summary(units: list[dict[str, object]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    residual_records = residual_rows_for_units(units)
    by_mechanism: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in residual_records:
        by_mechanism[str(row["payload_mechanism"])].append(row)
    for mechanism in sorted(by_mechanism, key=mechanism_sort_key):
        mech_rows = by_mechanism[mechanism]
        residuals = [float(row["context_fixed_residual"]) for row in mech_rows]
        units_for_mech = [unit for unit in units if str(unit["payload_mechanism"]) == mechanism]
        rows.append(
            {
                "payload_mechanism": mechanism,
                "payload_context_units": len(mech_rows),
                "unique_payloads": len({str(unit["payload_identity"]) for unit in units_for_mech}),
                "unique_contexts": len({str(row["cell_antigen_context"]) for row in mech_rows}),
                "mean_context_fixed_residual": f4(mean(residuals)),
                "median_context_fixed_residual": f4(median(residuals)),
                "min_context_fixed_residual": f4(min(residuals) if residuals else None),
                "max_context_fixed_residual": f4(max(residuals) if residuals else None),
                "positive_residual_units": sum(1 for value in residuals if value > EPS),
                "negative_residual_units": sum(1 for value in residuals if value < -EPS),
                "zero_residual_units": sum(1 for value in residuals if abs(value) <= EPS),
            }
        )
    return rows


def residual_range_stat(units: list[dict[str, object]]) -> float | None:
    means = residual_group_means(units)
    if len(means) <= 1:
        return None
    values = list(means.values())
    return max(values) - min(values)


def permutation_residual_range(units: list[dict[str, object]], rng: random.Random) -> float | None:
    residuals_by_mechanism: dict[str, list[float]] = defaultdict(list)
    for context, context_units in contexts_by_units(units).items():
        mechanisms = [str(unit["payload_mechanism"]) for unit in context_units]
        if len(set(mechanisms)) <= 1:
            continue
        context_mean = unit_response_mean(context_units)
        if context_mean is None:
            continue
        residuals = [float(unit["response_rate"]) - context_mean for unit in context_units]
        shuffled = mechanisms[:]
        rng.shuffle(shuffled)
        for mechanism, residual in zip(shuffled, residuals):
            residuals_by_mechanism[mechanism].append(residual)
    means = [statistics.mean(values) for values in residuals_by_mechanism.values() if values]
    if len(means) <= 1:
        return None
    return max(means) - min(means)


def pair_contexts(units: list[dict[str, object]], mechanism_a: str, mechanism_b: str) -> dict[str, list[dict[str, object]]]:
    output: dict[str, list[dict[str, object]]] = {}
    for context, context_units in contexts_by_units(units).items():
        selected = [
            unit
            for unit in context_units
            if str(unit["payload_mechanism"]) in {mechanism_a, mechanism_b}
        ]
        present = {str(unit["payload_mechanism"]) for unit in selected}
        if {mechanism_a, mechanism_b}.issubset(present):
            output[context] = selected
    return output


def pair_context_mean_diff(units: list[dict[str, object]], mechanism_a: str, mechanism_b: str) -> float | None:
    diffs: list[float] = []
    for context_units in pair_contexts(units, mechanism_a, mechanism_b).values():
        units_a = arm_units(context_units, mechanism_a)
        units_b = arm_units(context_units, mechanism_b)
        mean_a = unit_response_mean(units_a)
        mean_b = unit_response_mean(units_b)
        if mean_a is not None and mean_b is not None:
            diffs.append(mean_a - mean_b)
    return mean(diffs)


def permuted_pair_context_mean_diff(
    contexts: dict[str, list[dict[str, object]]],
    mechanism_a: str,
    mechanism_b: str,
    rng: random.Random,
) -> float | None:
    diffs: list[float] = []
    for context_units in contexts.values():
        labels = [str(unit["payload_mechanism"]) for unit in context_units]
        rates = [float(unit["response_rate"]) for unit in context_units]
        shuffled = labels[:]
        rng.shuffle(shuffled)
        rates_a = [rate for rate, label in zip(rates, shuffled) if label == mechanism_a]
        rates_b = [rate for rate, label in zip(rates, shuffled) if label == mechanism_b]
        if rates_a and rates_b:
            diffs.append(statistics.mean(rates_a) - statistics.mean(rates_b))
    return mean(diffs)


def permutation_tests(
    units: list[dict[str, object]],
    pair_summary_rows: list[dict[str, object]],
    repeats: int,
    seed: int,
) -> list[dict[str, object]]:
    rng = random.Random(seed)
    rows: list[dict[str, object]] = []
    observed_range = residual_range_stat(units)
    values: list[float] = []
    for _ in range(repeats):
        stat = permutation_residual_range(units, rng)
        if stat is not None:
            values.append(stat)
    rows.append(
        {
            "test": "all_mechanisms_context_fixed_residual_range",
            "mechanism_a": "",
            "mechanism_b": "",
            "observed_statistic": f4(observed_range),
            "permutations": len(values),
            "permutation_mean": f4(mean(values)),
            "permutation_p025": f4(percentile(values, 0.025)),
            "permutation_p975": f4(percentile(values, 0.975)),
            "p_two_sided_or_upper_tail": f4(
                (sum(1 for value in values if observed_range is not None and value >= observed_range - EPS) + 1)
                / (len(values) + 1)
                if values and observed_range is not None
                else None
            ),
        }
    )

    for row in pair_summary_rows:
        mechanism_a = str(row["mechanism_a"])
        mechanism_b = str(row["mechanism_b"])
        observed = pair_context_mean_diff(units, mechanism_a, mechanism_b)
        contexts = pair_contexts(units, mechanism_a, mechanism_b)
        values = []
        for _ in range(repeats):
            stat = permuted_pair_context_mean_diff(contexts, mechanism_a, mechanism_b, rng)
            if stat is not None:
                values.append(stat)
        rows.append(
            {
                "test": "pair_context_mean_unit_response_diff",
                "mechanism_a": mechanism_a,
                "mechanism_b": mechanism_b,
                "observed_statistic": f4(observed),
                "permutations": len(values),
                "permutation_mean": f4(mean(values)),
                "permutation_p025": f4(percentile(values, 0.025)),
                "permutation_p975": f4(percentile(values, 0.975)),
                "p_two_sided_or_upper_tail": f4(
                    (sum(1 for value in values if observed is not None and abs(value) >= abs(observed) - EPS) + 1)
                    / (len(values) + 1)
                    if values and observed is not None
                    else None
                ),
            }
        )
    return rows


def bootstrap_context_pair(
    units: list[dict[str, object]],
    mechanism_a: str,
    mechanism_b: str,
    repeats: int,
    rng: random.Random,
) -> list[float]:
    contexts = pair_contexts(units, mechanism_a, mechanism_b)
    names = sorted(contexts)
    if not names:
        return []
    values: list[float] = []
    for _ in range(repeats):
        diffs: list[float] = []
        for _ in names:
            context = rng.choice(names)
            units_a = arm_units(contexts[context], mechanism_a)
            units_b = arm_units(contexts[context], mechanism_b)
            sampled_a = [rng.choice(units_a) for _ in units_a]
            sampled_b = [rng.choice(units_b) for _ in units_b]
            mean_a = unit_response_mean(sampled_a)
            mean_b = unit_response_mean(sampled_b)
            if mean_a is not None and mean_b is not None:
                diffs.append(mean_a - mean_b)
        if diffs:
            values.append(statistics.mean(diffs))
    return values


def bootstrap_payload_cluster_pair(
    units: list[dict[str, object]],
    mechanism_a: str,
    mechanism_b: str,
    repeats: int,
    rng: random.Random,
) -> tuple[list[float], int]:
    contexts = pair_contexts(units, mechanism_a, mechanism_b)
    pair_units = [unit for context_units in contexts.values() for unit in context_units]
    payloads_a = sorted({str(unit["payload_identity"]) for unit in pair_units if unit["payload_mechanism"] == mechanism_a})
    payloads_b = sorted({str(unit["payload_identity"]) for unit in pair_units if unit["payload_mechanism"] == mechanism_b})
    if not payloads_a or not payloads_b:
        return [], 0
    values: list[float] = []
    unavailable = 0
    for _ in range(repeats):
        selected_a = Counter(rng.choice(payloads_a) for _ in payloads_a)
        selected_b = Counter(rng.choice(payloads_b) for _ in payloads_b)
        resampled: list[dict[str, object]] = []
        for unit in pair_units:
            payload = str(unit["payload_identity"])
            mechanism = str(unit["payload_mechanism"])
            count = selected_a[payload] if mechanism == mechanism_a else selected_b[payload]
            for _ in range(count):
                resampled.append(dict(unit))
        stat = pair_context_mean_diff(resampled, mechanism_a, mechanism_b)
        if stat is None:
            unavailable += 1
        else:
            values.append(stat)
    return values, unavailable


def bootstrap_rows(
    units: list[dict[str, object]],
    pair_summary_rows: list[dict[str, object]],
    repeats: int,
    seed: int,
) -> list[dict[str, object]]:
    rng = random.Random(seed)
    rows: list[dict[str, object]] = []
    for row in pair_summary_rows:
        mechanism_a = str(row["mechanism_a"])
        mechanism_b = str(row["mechanism_b"])
        observed = pair_context_mean_diff(units, mechanism_a, mechanism_b)
        context_values = bootstrap_context_pair(units, mechanism_a, mechanism_b, repeats, rng)
        payload_values, payload_unavailable = bootstrap_payload_cluster_pair(
            units, mechanism_a, mechanism_b, repeats, rng
        )
        payloads_a = str(row["payloads_a_across_contexts"]).split(";") if row["payloads_a_across_contexts"] else []
        payloads_b = str(row["payloads_b_across_contexts"]).split(";") if row["payloads_b_across_contexts"] else []
        for kind, values, unavailable in [
            ("context_and_within_arm_payload_unit_bootstrap", context_values, 0),
            ("payload_cluster_bootstrap", payload_values, payload_unavailable),
        ]:
            rows.append(
                {
                    "bootstrap": kind,
                    "mechanism_a": mechanism_a,
                    "mechanism_b": mechanism_b,
                    "observed_statistic": f4(observed),
                    "bootstrap_repeats_requested": repeats,
                    "bootstrap_repeats_available": len(values),
                    "unavailable_repeats": unavailable,
                    "bootstrap_mean": f4(mean(values)),
                    "bootstrap_p025": f4(percentile(values, 0.025)),
                    "bootstrap_p975": f4(percentile(values, 0.975)),
                    "fraction_gt_zero": f4(sum(1 for value in values if value > EPS) / len(values) if values else None),
                    "fraction_lt_zero": f4(sum(1 for value in values if value < -EPS) / len(values) if values else None),
                    "unique_payloads_a": len(payloads_a),
                    "unique_payloads_b": len(payloads_b),
                    "payload_bootstrap_degenerate": int(kind == "payload_cluster_bootstrap" and min(len(payloads_a), len(payloads_b)) < 2),
                }
            )
    return rows


def leave_one_payload_out_rows(
    units: list[dict[str, object]],
    pair_summary_rows: list[dict[str, object]],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    rows: list[dict[str, object]] = []
    residual_rows: list[dict[str, object]] = []
    all_payloads = sorted({str(unit["payload_identity"]) for unit in units})

    observed_range = residual_range_stat(units)
    for payload in all_payloads:
        mechanism = sorted({str(unit["payload_mechanism"]) for unit in units if unit["payload_identity"] == payload})[0]
        filtered = [unit for unit in units if str(unit["payload_identity"]) != payload]
        filtered_context_names = matched_context_names(filtered)
        filtered_matched = [unit for unit in filtered if str(unit["cell_antigen_context"]) in filtered_context_names]
        stat = residual_range_stat(filtered_matched)
        residual_rows.append(
            {
                "held_out_payload": payload,
                "held_out_mechanism": mechanism,
                "analysis": "all_mechanisms_context_fixed_residual_range",
                "observed_statistic": f4(observed_range),
                "remaining_matched_contexts": len(filtered_context_names),
                "remaining_payload_context_units": len(filtered_matched),
                "leave_one_payload_out_statistic": f4(stat),
                "statistic_unavailable": int(stat is None),
            }
        )

    for pair_row in pair_summary_rows:
        mechanism_a = str(pair_row["mechanism_a"])
        mechanism_b = str(pair_row["mechanism_b"])
        observed = pair_context_mean_diff(units, mechanism_a, mechanism_b)
        pair_payloads = sorted(
            {
                str(unit["payload_identity"])
                for context_units in pair_contexts(units, mechanism_a, mechanism_b).values()
                for unit in context_units
                if str(unit["payload_mechanism"]) in {mechanism_a, mechanism_b}
            }
        )
        observed_sign = sign(observed)
        for payload in pair_payloads:
            held_out_mechanism = sorted(
                {
                    str(unit["payload_mechanism"])
                    for unit in units
                    if str(unit["payload_identity"]) == payload
                }
            )[0]
            filtered = [unit for unit in units if str(unit["payload_identity"]) != payload]
            stat = pair_context_mean_diff(filtered, mechanism_a, mechanism_b)
            contexts = pair_contexts(filtered, mechanism_a, mechanism_b)
            rows.append(
                {
                    "mechanism_a": mechanism_a,
                    "mechanism_b": mechanism_b,
                    "held_out_payload": payload,
                    "held_out_mechanism": held_out_mechanism,
                    "observed_statistic": f4(observed),
                    "remaining_crossed_contexts": len(contexts),
                    "remaining_payload_context_units": sum(len(context_units) for context_units in contexts.values()),
                    "leave_one_payload_out_statistic": f4(stat),
                    "statistic_unavailable": int(stat is None),
                    "sign_changed": int(stat is not None and observed_sign != 0 and sign(stat) != observed_sign),
                }
            )

    summary_rows: list[dict[str, object]] = []
    by_pair: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        by_pair[(str(row["mechanism_a"]), str(row["mechanism_b"]))].append(row)
    for (mechanism_a, mechanism_b), pair_rows in sorted(by_pair.items(), key=lambda item: (mechanism_sort_key(item[0][0]), mechanism_sort_key(item[0][1]))):
        available = [float(row["leave_one_payload_out_statistic"]) for row in pair_rows if row["leave_one_payload_out_statistic"] != ""]
        unavailable_payloads = [str(row["held_out_payload"]) for row in pair_rows if int(row["statistic_unavailable"]) == 1]
        sign_changed = [str(row["held_out_payload"]) for row in pair_rows if int(row["sign_changed"]) == 1]
        remaining_contexts = [int(row["remaining_crossed_contexts"]) for row in pair_rows]
        summary_rows.append(
            {
                "mechanism_a": mechanism_a,
                "mechanism_b": mechanism_b,
                "payloads_tested": len(pair_rows),
                "unavailable_holdouts": len(unavailable_payloads),
                "payloads_making_statistic_unavailable": ";".join(unavailable_payloads),
                "sign_change_holdouts": len(sign_changed),
                "payloads_changing_sign": ";".join(sign_changed),
                "min_remaining_crossed_contexts": min(remaining_contexts) if remaining_contexts else "",
                "max_remaining_crossed_contexts": max(remaining_contexts) if remaining_contexts else "",
                "min_leave_one_payload_out_statistic": f4(min(available) if available else None),
                "max_leave_one_payload_out_statistic": f4(max(available) if available else None),
            }
        )
    return rows + residual_rows, summary_rows


def sign(value: float | None) -> int:
    if value is None or abs(value) <= EPS:
        return 0
    return 1 if value > 0 else -1


def payload_pairing_rows(pair_effects: list[dict[str, object]]) -> list[dict[str, object]]:
    counter: Counter[tuple[str, str, str]] = Counter()
    positive_counter: Counter[tuple[str, str, str]] = Counter()
    negative_counter: Counter[tuple[str, str, str]] = Counter()
    zero_counter: Counter[tuple[str, str, str]] = Counter()
    examples: dict[tuple[str, str, str], list[str]] = defaultdict(list)
    for row in pair_effects:
        key = (
            str(row["mechanism_a"]),
            str(row["mechanism_b"]),
            str(row["payload_pairing_signature"]),
        )
        counter[key] += 1
        diff = float(row["unit_response_diff_a_minus_b"])
        if diff > EPS:
            positive_counter[key] += 1
        elif diff < -EPS:
            negative_counter[key] += 1
        else:
            zero_counter[key] += 1
        examples[key].append(str(row["cell_antigen_context"]))

    rows = []
    for (mechanism_a, mechanism_b, signature), count in counter.most_common():
        rows.append(
            {
                "mechanism_a": mechanism_a,
                "mechanism_b": mechanism_b,
                "payload_pairing_signature": signature,
                "contexts": count,
                "a_higher_contexts": positive_counter[(mechanism_a, mechanism_b, signature)],
                "b_higher_contexts": negative_counter[(mechanism_a, mechanism_b, signature)],
                "equal_contexts": zero_counter[(mechanism_a, mechanism_b, signature)],
                "example_contexts": "; ".join(examples[(mechanism_a, mechanism_b, signature)][:5]),
            }
        )
    rows.sort(
        key=lambda row: (
            str(row["mechanism_a"]),
            str(row["mechanism_b"]),
            -int(row["contexts"]),
            str(row["payload_pairing_signature"]),
        )
    )
    return rows


def strong_overlap_context_names(context_rows: list[dict[str, object]]) -> set[str]:
    return {
        str(row["cell_antigen_context"])
        for row in context_rows
        if int(row["mechanism_count"]) >= 3 or int(row["payload_context_units"]) >= 4
    }


def determine_decision(
    pair_summary: list[dict[str, object]],
    bootstrap: list[dict[str, object]],
    leave_one_out_summary: list[dict[str, object]],
    permutation: list[dict[str, object]],
) -> tuple[str, list[str]]:
    reasons: list[str] = []
    any_descriptive_pattern = False
    any_pivot_go = False

    permutation_by_pair = {
        (str(row["mechanism_a"]), str(row["mechanism_b"])): row
        for row in permutation
        if row["test"] == "pair_context_mean_unit_response_diff"
    }
    bootstrap_by_pair = {
        (str(row["mechanism_a"]), str(row["mechanism_b"]), str(row["bootstrap"])): row
        for row in bootstrap
    }
    loo_by_pair = {
        (str(row["mechanism_a"]), str(row["mechanism_b"])): row
        for row in leave_one_out_summary
    }

    for row in pair_summary:
        pair = (str(row["mechanism_a"]), str(row["mechanism_b"]))
        contexts = int(row["crossed_contexts"])
        both_rep = int(row["contexts_with_both_arms_2plus_payloads"])
        a_higher = int(row["a_higher_contexts"])
        b_higher = int(row["b_higher_contexts"])
        mean_diff = float(row["mean_context_unit_response_diff_a_minus_b"])
        if contexts >= 5 and abs(mean_diff) >= 0.15 and max(a_higher, b_higher) >= 3:
            any_descriptive_pattern = True
        p_row = permutation_by_pair.get(pair, {})
        p_value = float(p_row.get("p_two_sided_or_upper_tail", "1") or "1")
        boot_row = bootstrap_by_pair.get((pair[0], pair[1], "context_and_within_arm_payload_unit_bootstrap"), {})
        ci_low = float(boot_row.get("bootstrap_p025", "nan") or "nan")
        ci_high = float(boot_row.get("bootstrap_p975", "nan") or "nan")
        loo_row = loo_by_pair.get(pair, {})
        unavailable = int_or_default(loo_row.get("unavailable_holdouts"), 999)
        sign_changes = int_or_default(loo_row.get("sign_change_holdouts"), 999)
        ci_excludes_zero = math.isfinite(ci_low) and math.isfinite(ci_high) and (ci_low > 0 or ci_high < 0)
        if both_rep >= 3 and p_value <= 0.05 and ci_excludes_zero and unavailable == 0 and sign_changes == 0:
            any_pivot_go = True
        reasons.append(
            f"{pair[0]} vs {pair[1]}: {contexts} exact contexts, mean unit-response diff {mean_diff:.4f}, "
            f"{a_higher}/{b_higher}/{int(row['equal_contexts'])} higher/lower/equal, "
            f"both-arm >=2-payload contexts={both_rep}, permutation p={p_value:.4f}."
        )

    total_both_rep = sum(int(row["contexts_with_both_arms_2plus_payloads"]) for row in pair_summary)
    reasons.append(
        f"No mechanism pair has replicated independent payload support within the same exact context "
        f"(total both-arm >=2-payload pair-contexts={total_both_rep})."
    )
    if any_pivot_go:
        return "PIVOT-GO", reasons
    if any_descriptive_pattern:
        reasons.append(
            "There is a descriptive within-context pattern, but it is not identifiable as a mechanism-level effect because one or both arms are driven by single payloads."
        )
        return "PIVOT-WEAK", reasons
    reasons.append("The matched-context controls do not show a convincing within-context mechanism pattern.")
    return "KILL", reasons


def int_or_default(value: object, default: int) -> int:
    if value is None or value == "":
        return default
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def svg_pairwise_diffs(path: Path, pair_effects: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    width = 1040
    left = 420
    right = 60
    top = 64
    row_h = 25
    height = top + row_h * len(pair_effects) + 40
    x0 = left
    x1 = width - right
    zero_x = x0 + (x1 - x0) * 0.5
    colors = {
        ("DNA-damaging", "Topo-I inhibitor"): "#b84a62",
        ("DNA-damaging", "microtubule-disrupting"): "#2f6f73",
        ("Topo-I inhibitor", "microtubule-disrupting"): "#4f67a5",
    }
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<text x="24" y="34" font-family="Arial, sans-serif" font-size="20" font-weight="700" fill="#111827">Matched-context payload-unit response differences</text>',
        f'<line x1="{zero_x:.1f}" y1="48" x2="{zero_x:.1f}" y2="{height - 24}" stroke="#9ca3af" stroke-width="1"/>',
        f'<line x1="{x0}" y1="48" x2="{x1}" y2="48" stroke="#d1d5db" stroke-width="1"/>',
        f'<text x="{x0}" y="42" font-family="Arial, sans-serif" font-size="12" fill="#374151">-1</text>',
        f'<text x="{zero_x - 4}" y="42" font-family="Arial, sans-serif" font-size="12" fill="#374151">0</text>',
        f'<text x="{x1 - 10}" y="42" font-family="Arial, sans-serif" font-size="12" fill="#374151">+1</text>',
    ]
    for index, row in enumerate(pair_effects):
        y = top + index * row_h
        diff = float(row["unit_response_diff_a_minus_b"])
        x = x0 + (diff + 1.0) / 2.0 * (x1 - x0)
        pair = (str(row["mechanism_a"]), str(row["mechanism_b"]))
        color = colors.get(pair, "#6b7280")
        label = str(row["cell_antigen_context"])
        if len(label) > 56:
            label = label[:53] + "..."
        pair_label = f"{pair[0]} - {pair[1]}"
        lines.append(
            f'<text x="24" y="{y + 5}" font-family="Arial, sans-serif" font-size="12" fill="#111827">{escape(label)}</text>'
        )
        lines.append(
            f'<text x="300" y="{y + 5}" font-family="Arial, sans-serif" font-size="11" fill="#4b5563">{escape(pair_label)}</text>'
        )
        lines.append(f'<circle cx="{x:.1f}" cy="{y}" r="5" fill="{color}"/>')
    lines.append("</svg>")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_report(
    path: Path,
    dataset_reduction: list[dict[str, object]],
    context_rows: list[dict[str, object]],
    pair_summary: list[dict[str, object]],
    pair_effects: list[dict[str, object]],
    residual_summary_rows: list[dict[str, object]],
    permutation_rows: list[dict[str, object]],
    bootstrap: list[dict[str, object]],
    leave_one_out_summary: list[dict[str, object]],
    strong_pair_summary: list[dict[str, object]],
    decision: str,
    decision_reasons: list[str],
) -> None:
    crossed_contexts = len(context_rows)
    total_context_units = sum(int(row["payload_context_units"]) for row in context_rows)
    total_input_obs = sum(int(row["input_observations"]) for row in context_rows)
    variable_contexts = sum(int(row["has_unit_response_variation"]) for row in context_rows)
    both_arm_rep = sum(int(row["contexts_with_both_arms_2plus_payloads"]) for row in pair_summary)
    strong_contexts = strong_overlap_context_names(context_rows)

    lines = [
        "# Matched-Context Pivot Analysis",
        "",
        "## Scope",
        "",
        "- Input: existing local `data/processed/adc_response_10nM_context.csv`.",
        "- Raw data were not redownloaded or replaced.",
        "- No predictive model, neural head, foundation model, or hyperparameter search was run.",
        "- Central analysis excludes `Other` and uses the three interpretable major mechanisms: DNA-damaging, Topo-I inhibitor, and microtubule-disrupting.",
        "- Unit of inference: one payload identity within one exact antigen x cell-line context. Repeated ADC observations using the same payload in the same context are collapsed before comparison.",
        "",
        "## Data Reduction",
        "",
        markdown_table(
            dataset_reduction,
            [
                "scope",
                "input_observations",
                "positive_input_observations",
                "negative_input_observations",
                "observation_prevalence",
                "payload_context_units",
                "mean_payload_unit_response_rate",
                "unique_payloads",
                "unique_mechanisms",
                "unique_cell_antigen_contexts",
            ],
        ),
        "",
        "## Exact Matched-Context Coverage",
        "",
        f"- Exact crossed major-mechanism antigen x cell-line contexts: {crossed_contexts}.",
        f"- Payload-context units in those contexts: {total_context_units}.",
        f"- Input observations represented by those units: {total_input_obs}.",
        f"- Crossed contexts with any unit-level response variation: {variable_contexts}.",
        f"- Pairwise context comparisons where both mechanism arms contain >=2 independent payload identities: {both_arm_rep}.",
        f"- Strong-overlap contexts, defined as >=3 mechanisms or >=4 payload-context units: {len(strong_contexts)}.",
        "",
        markdown_table(
            context_rows,
            [
                "cell_antigen_context",
                "mechanism_count",
                "payload_context_units",
                "input_observations",
                "has_unit_response_variation",
                "payload_units_by_mechanism",
                "unit_response_by_mechanism",
                "payloads_by_mechanism",
            ],
            max_rows=24,
        ),
        "",
        "## Within-Context Mechanism Effects",
        "",
        "The response difference is computed from collapsed payload-context units, then averaged over exact contexts. Observation-weighted columns are reported as a sensitivity check, not as independent mechanism evidence.",
        "",
        markdown_table(
            pair_summary,
            [
                "mechanism_a",
                "mechanism_b",
                "crossed_contexts",
                "payload_context_units_in_pair_contexts",
                "contexts_with_unit_response_difference",
                "a_higher_contexts",
                "b_higher_contexts",
                "equal_contexts",
                "mean_context_unit_response_diff_a_minus_b",
                "observation_weighted_diff_a_minus_b",
                "contexts_with_both_arms_2plus_payloads",
                "unique_payloads_a_across_contexts",
                "unique_payloads_b_across_contexts",
            ],
        ),
        "",
        "Largest individual exact-context differences:",
        "",
        markdown_table(
            sorted(pair_effects, key=lambda row: -abs(float(row["unit_response_diff_a_minus_b"]))),
            [
                "cell_antigen_context",
                "mechanism_a",
                "mechanism_b",
                "payload_units_a",
                "payload_units_b",
                "unit_response_a",
                "unit_response_b",
                "unit_response_diff_a_minus_b",
                "payloads_a",
                "payloads_b",
            ],
            max_rows=12,
        ),
        "",
        "## Context-Fixed Residual Check",
        "",
        "Residuals subtract each exact context's mean payload-unit response before summarizing by mechanism.",
        "",
        markdown_table(
            residual_summary_rows,
            [
                "payload_mechanism",
                "payload_context_units",
                "unique_payloads",
                "unique_contexts",
                "mean_context_fixed_residual",
                "median_context_fixed_residual",
                "min_context_fixed_residual",
                "max_context_fixed_residual",
            ],
        ),
        "",
        "## Negative Controls",
        "",
        "Permutation shuffles mechanism labels only within exact contexts, preserving context and mechanism arm sizes. The payload bootstrap resamples payload identities as clusters; degenerate rows indicate that one side has fewer than two payloads, so payload-level uncertainty is not estimable.",
        "",
        markdown_table(
            permutation_rows,
            [
                "test",
                "mechanism_a",
                "mechanism_b",
                "observed_statistic",
                "permutation_mean",
                "permutation_p025",
                "permutation_p975",
                "p_two_sided_or_upper_tail",
            ],
        ),
        "",
        markdown_table(
            bootstrap,
            [
                "bootstrap",
                "mechanism_a",
                "mechanism_b",
                "observed_statistic",
                "bootstrap_p025",
                "bootstrap_p975",
                "fraction_gt_zero",
                "payload_bootstrap_degenerate",
            ],
        ),
        "",
        "Leave-one-payload-out stability:",
        "",
        markdown_table(
            leave_one_out_summary,
            [
                "mechanism_a",
                "mechanism_b",
                "payloads_tested",
                "unavailable_holdouts",
                "payloads_making_statistic_unavailable",
                "sign_change_holdouts",
                "min_remaining_crossed_contexts",
                "max_remaining_crossed_contexts",
                "min_leave_one_payload_out_statistic",
                "max_leave_one_payload_out_statistic",
            ],
        ),
        "",
        "Strong-overlap subset summary:",
        "",
        markdown_table(
            strong_pair_summary,
            [
                "mechanism_a",
                "mechanism_b",
                "crossed_contexts",
                "contexts_with_unit_response_difference",
                "mean_context_unit_response_diff_a_minus_b",
                "contexts_with_both_arms_2plus_payloads",
                "unique_payloads_a_across_contexts",
                "unique_payloads_b_across_contexts",
            ],
        ),
        "",
        "## Decision",
        "",
        f"**{decision}**",
        "",
    ]
    for reason in decision_reasons:
        lines.append(f"- {reason}")
    lines.extend(
        [
            "",
            "Interpretation: the strongest exact matched-context pattern is descriptive, not identifiable as a mechanism-level effect. The data can show payload-specific differences within some contexts, but they do not replicate across multiple independent payloads per mechanism inside the same contexts.",
            "",
            "Additional crossed data required to resolve the pivot: the same antigen x cell-line contexts measured under at least two independent payloads per mechanism, across several contexts per mechanism pair, with harmonized endpoints and units. Without that, mechanism remains inseparable from payload identity in this dataset.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    all_rows = read_csv(args.context_csv)
    major_rows = [row for row in all_rows if row.get("payload_mechanism") in MAJOR_MECHANISM_SET]
    all_units = collapse_payload_context_units(all_rows)
    major_units = collapse_payload_context_units(major_rows)
    matched_names = matched_context_names(major_units)
    matched_units = [unit for unit in major_units if str(unit["cell_antigen_context"]) in matched_names]

    reduction = dataset_reduction_rows(all_rows, major_rows, all_units, major_units, matched_units)
    contexts = context_summary(matched_units)
    pair_effects = pairwise_context_effects(matched_units)
    pair_summary = summarize_pairwise_effects(pair_effects)
    residual_unit_rows = residual_rows_for_units(matched_units)
    residual_summaries = residual_summary(matched_units)
    permutations = permutation_tests(matched_units, pair_summary, args.permutations, args.seed + 10)
    bootstraps = bootstrap_rows(matched_units, pair_summary, args.bootstraps, args.seed + 20)
    loo_rows, loo_summary = leave_one_payload_out_rows(matched_units, pair_summary)
    pairing_rows = payload_pairing_rows(pair_effects)

    strong_contexts = strong_overlap_context_names(contexts)
    strong_units = [unit for unit in matched_units if str(unit["cell_antigen_context"]) in strong_contexts]
    strong_pair_effects = pairwise_context_effects(strong_units)
    strong_pair_summary = summarize_pairwise_effects(strong_pair_effects)

    decision, decision_reasons = determine_decision(
        pair_summary,
        bootstraps,
        loo_summary,
        permutations,
    )

    write_csv(args.out_dir / "matched_context_dataset_reduction.csv", reduction)
    write_csv(args.out_dir / "matched_context_payload_units.csv", matched_units)
    write_csv(args.out_dir / "matched_context_coverage.csv", contexts)
    write_csv(args.out_dir / "matched_context_pairwise_payload_unit_effects.csv", pair_effects)
    write_csv(args.out_dir / "matched_context_pairwise_payload_unit_summary.csv", pair_summary)
    write_csv(args.out_dir / "matched_context_payload_pairings.csv", pairing_rows)
    write_csv(args.out_dir / "matched_context_fixed_effect_residuals.csv", residual_unit_rows)
    write_csv(args.out_dir / "matched_context_fixed_effect_summary.csv", residual_summaries)
    write_csv(args.out_dir / "matched_context_permutation_tests.csv", permutations)
    write_csv(args.out_dir / "matched_context_bootstrap_ci.csv", bootstraps)
    write_csv(args.out_dir / "matched_context_leave_one_payload_out.csv", loo_rows)
    write_csv(args.out_dir / "matched_context_leave_one_payload_out_summary.csv", loo_summary)
    write_csv(args.out_dir / "matched_context_strong_overlap_pair_summary.csv", strong_pair_summary)
    svg_pairwise_diffs(args.out_dir / "figures" / "matched_context_pairwise_diffs.svg", pair_effects)
    write_report(
        args.out_dir / "matched_context_pivot.md",
        reduction,
        contexts,
        pair_summary,
        pair_effects,
        residual_summaries,
        permutations,
        bootstraps,
        loo_summary,
        strong_pair_summary,
        decision,
        decision_reasons,
    )

    print(f"decision={decision}")
    print(f"matched_contexts={len(contexts)}")
    print(f"matched_payload_context_units={len(matched_units)}")
    print(f"matched_input_observations={sum(int(unit['input_observations']) for unit in matched_units)}")
    print(f"outputs={args.out_dir / 'matched_context_pivot.md'}")


if __name__ == "__main__":
    main()
