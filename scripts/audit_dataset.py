#!/usr/bin/env python3
"""Minimal ADC dataset audit for the mechanism-aware pilot."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


DEFAULT_COLUMNS = {
    "mechanism": "payload_mechanism",
    "payload": "payload",
    "antigen": "antigen",
    "cell_line": "cell_line",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv", type=Path, help="Input CSV")
    parser.add_argument("--mechanism-col", default=DEFAULT_COLUMNS["mechanism"])
    parser.add_argument("--payload-col", default=DEFAULT_COLUMNS["payload"])
    parser.add_argument("--antigen-col", default=DEFAULT_COLUMNS["antigen"])
    parser.add_argument("--cell-line-col", default=DEFAULT_COLUMNS["cell_line"])
    parser.add_argument("--out", type=Path, default=Path("results/data_audit.md"))
    return parser.parse_args()


def safe_nunique(df: pd.DataFrame, col: str) -> int | None:
    return int(df[col].nunique(dropna=True)) if col in df.columns else None


def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.csv)

    if args.mechanism_col not in df.columns:
        raise ValueError(f"Missing required mechanism column: {args.mechanism_col}")

    mech = args.mechanism_col
    counts = df[mech].fillna("<missing>").value_counts(dropna=False)

    lines = [
        "# ADC dataset audit",
        "",
        f"- Rows: **{len(df):,}**",
        f"- Mechanism classes: **{df[mech].nunique(dropna=True)}**",
        f"- Unique payloads: **{safe_nunique(df, args.payload_col)}**",
        f"- Unique antigens: **{safe_nunique(df, args.antigen_col)}**",
        f"- Unique cell lines: **{safe_nunique(df, args.cell_line_col)}**",
        "",
        "## Observations per mechanism",
        "",
        counts.to_frame("n").to_markdown(),
        "",
    ]

    if args.payload_col in df.columns:
        payload_by_mech = df.groupby(mech, dropna=False)[args.payload_col].nunique(dropna=True)
        lines += ["## Unique payloads per mechanism", "", payload_by_mech.to_frame("unique_payloads").to_markdown(), ""]

        pairs = df[[mech, args.payload_col]].dropna().drop_duplicates()
        payload_mech_degree = pairs.groupby(args.payload_col)[mech].nunique()
        one_to_one_fraction = float((payload_mech_degree == 1).mean()) if len(payload_mech_degree) else float("nan")
        lines += [
            "## Mechanism–payload confounding diagnostic",
            "",
            f"Fraction of payloads observed in exactly one mechanism class: **{one_to_one_fraction:.3f}**",
            "",
            "A value near 1.0 indicates strong mechanism/payload confounding and limits claims that a model learned transferable mechanism biology rather than payload identity.",
            "",
        ]

    for label, col in [("antigens", args.antigen_col), ("cell lines", args.cell_line_col)]:
        if col in df.columns:
            overlap = df.groupby(col)[mech].nunique(dropna=True)
            multi = int((overlap > 1).sum())
            lines += [
                f"## Cross-mechanism overlap: {label}",
                "",
                f"{multi:,} / {len(overlap):,} unique {label} occur in more than one mechanism class.",
                "",
            ]

    lines += [
        "## Decision",
        "",
        "- GO / KILL / PIVOT: **TBD**",
        "- Main identifiability concern: **TBD**",
        "- Cheapest next discriminating experiment: **TBD**",
        "",
    ]

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
