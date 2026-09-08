# ADCMechanism

Mechanism-aware antibody–drug conjugate (ADC) response modeling.

## Research question

Can explicit modeling of payload mechanism improve ADC response prediction and biological interpretability by separating target/delivery biology from payload-specific sensitivity?

## Initial falsification pilot

Before training a complex model, audit whether the available data can support the hypothesis.

1. Count observations per payload-mechanism class.
2. Count unique payloads, antigens, and cell lines per mechanism.
3. Measure class imbalance and mechanism–payload confounding.
4. Measure overlap of biological contexts across mechanisms.
5. Decide whether a mechanism-aware model is identifiable from the available data.

Only proceed to modeling if the audit passes.

## Modeling ladder

1. Simple response baseline.
2. Baseline + mechanism indicator.
3. Shared encoder + mechanism-specific heads.
4. Explicit delivery / payload-sensitivity feature separation.
5. OOD evaluation, prioritizing held-out payloads or mechanism-relevant generalization.

## Repository layout

```text
ADCMechanism/
├── configs/        # experiment configuration
├── data/           # local/raw/intermediate data (not committed)
├── docs/           # notes and decisions
├── results/        # compact tables/figures/metrics
├── scripts/        # executable analyses and experiments
├── src/            # reusable project code
├── tests/          # lightweight checks
├── .gitignore
├── pyproject.toml
└── README.md
```

## Research operating rule

Prefer the cheapest experiment that can falsify the current hypothesis. Do not build the full model or benchmark before the data audit supports doing so.

For each meaningful run, record:

- what was tested
- what was observed
- whether the hypothesis became stronger or weaker
- the next cheapest discriminating experiment
- **GO / KILL / PIVOT**

## Getting started

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .
python scripts/audit_dataset.py --help
```

## Current local study

The current ADCdb/ADCNet mechanism identifiability run is reproducible with:

```bash
python3 scripts/audit_adcdb_activity.py
python3 scripts/identifiability_study.py
python3 scripts/matched_context_pivot.py
```

Primary outputs:

- `results/identifiability_study.md`
- `results/cheap_model_results.csv`
- `results/matched_context_mechanism_effect_summary.csv`
- `results/matched_context_pivot.md`
- `OVERNIGHT_STATUS.md`
