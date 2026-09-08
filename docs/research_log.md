# Research Log

Use one dated section per meaningful experiment or decision.

## Template

### YYYY-MM-DD — Short experiment name

**Question**

What specific uncertainty is this experiment trying to reduce?

**Experiment**

What was run, with enough detail to reproduce it?

**Result**

What happened? Include the key quantitative result or failure mode.

**Interpretation**

Did this strengthen or weaken the mechanism-aware ADC hypothesis?

**Decision**

**GO / KILL / PIVOT**

**Next cheapest discriminating experiment**

What is the smallest next test that could change the decision?

### 2026-09-08 - ADCdb mechanism identifiability and cheap falsification

**Question**

Can the local ADCdb/ADCNet-derived response data support a defensible test that explicit payload mechanism improves response prediction beyond payload identity, antigen/target, and measured cell-line context?

**Experiment**

Commands:

```bash
python3 scripts/audit_adcdb_activity.py
python3 scripts/identifiability_study.py
```

The second command read local `data/raw/` files, regenerated compact processed tables, derived deterministic 10 nM labels from molar IC50/LC50/EC50/GI50 rows, collapsed duplicate ADC/cell-line/antigen/payload rows, audited mechanism confounding, and ran categorical Naive Bayes metadata baselines. No neural model was trained.

Primary modeling dataset:

- 1,043 parsed ADCdb cell-line activity rows.
- 705 usable deterministic molar potency rows at 10 nM.
- 624 collapsed observations; 437 positives and 187 negatives.
- Major-mechanism subset for modeling excludes `Other`: 609 observations; 428 positives and 181 negatives.

Splits:

- `random_row`: stratified row split; included only to expose shortcut behavior.
- `payload_group`: 5-fold grouped by payload identity; central split.
- `payload_leave_one_out`: one payload held out per fold.
- `cell_line_group`: grouped by cell line.
- `cell_antigen_context_group`: grouped by cell-line/antigen context.

**Result**

Phase 2 audit was sufficient for a cheap falsification run, but not for a full architecture:

- Major mechanisms each have at least 6 payloads and at least 78 observations.
- 24 crossed major-mechanism cell-line/antigen contexts cover 113 observations.
- Only 6 crossed cell-line/antigen contexts contain both response classes.
- Mechanism is perfectly associated with payload identity in the collapsed table: Cramer's V = 1.000.
- Mechanism is also strongly associated with antigen, cell line, and cell-line/antigen context: Cramer's V = 0.828, 0.877, and 0.942.

Central payload-held-out result:

- `bio` (antigen + cell line): AP 0.8341, ROC-AUC 0.6805, MCC 0.2788.
- `bio_plus_mechanism`: AP 0.8294, ROC-AUC 0.6819, MCC 0.2159.
- `bio_plus_payload`: AP 0.8422, ROC-AUC 0.6931, MCC 0.2562.
- `bio_plus_payload_plus_mechanism`: AP 0.8446, ROC-AUC 0.6949, MCC 0.1852.

Stress checks:

- Payload leave-one-out also did not support a mechanism gain: `bio` AP 0.8333 vs `bio_plus_mechanism` AP 0.8203.
- Payload-level shuffled mechanism control under `payload_group`: true `bio_plus_mechanism` AP 0.8294; shuffled mean AP 0.8261; shuffled p95 AP 0.8497.
- Matched cell-line/antigen context effects were sparse: 7/28 pairwise contexts had label variation.

**Interpretation**

The local data contain enough structure for a cheap falsification experiment, but the result does not support the original mechanism-aware hypothesis. Payload mechanism is not adding robust held-out-payload signal beyond antigen/cell-line metadata, and any apparent mechanism signal is comparable to shuffled mechanism labels. Strong mechanism associations with payload, antigen, cell line, and cell-line/antigen context make random-split gains scientifically weak.

**Decision**

**PIVOT**

The defensible reformulation is narrower: test whether coarse payload mechanism can act as a low-dimensional prior for unseen payloads in ADCdb-like metadata, not whether it independently explains ADC response beyond payload identity and biological context.

**Next cheapest discriminating experiment**

Stop model escalation. Curate before modeling: validate cell-line naming/synonyms, endpoint aggregation, and mass-unit conversion assumptions to see whether additional genuinely matched cell-line/antigen contexts can be recovered from the existing local ADCdb pages.
