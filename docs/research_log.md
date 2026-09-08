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

### 2026-09-08 - Exact matched-context pivot analysis

**Question**

After controlling biological context as tightly as the current data permit, is payload mechanism associated with differential ADC response within the same antigen x cell-line contexts?

**Experiment**

Command:

```bash
python3 scripts/matched_context_pivot.py
```

The script reads the existing processed ADCdb table only. It excludes `Other`, collapses repeated observations to one payload-context unit per payload identity within each exact antigen x cell-line context, identifies contexts measured under multiple major mechanisms, and runs transparent within-context comparisons, context-fixed residual summaries, within-context permutation controls, context/payload bootstraps, leave-one-payload-out checks, and a strongest-overlap subset analysis. No predictive model was trained.

**Result**

- 24 exact crossed major-mechanism antigen x cell-line contexts were available.
- These represented 113 input observations, but only 58 payload-context units after collapsing repeated payload observations.
- The matched subset retained 18.6% of major-mechanism input observations and 13.2% of major-mechanism payload-context units.
- 6/24 exact crossed contexts had any unit-level response variation.
- 0 pairwise context comparisons had >=2 independent payload identities in both mechanism arms.

Pairwise payload-unit response differences:

- DNA-damaging vs microtubule-disrupting: 16 contexts, mean difference 0.0030, 1/1/14 higher/lower/equal, permutation p=1.0000.
- DNA-damaging vs Topo-I inhibitor: 7 contexts, mean difference 0.3571, 4/0/3 higher/lower/equal, permutation p=0.1266.
- Topo-I inhibitor vs microtubule-disrupting: 5 contexts, mean difference -0.2000, 0/1/4 higher/lower/equal, permutation p=0.3294.

The strongest descriptive signal was DNA-damaging higher than Topo-I, but Topo-I was represented by only `PAY0ZVBAI` in that comparison. Leave-one-payload-out confirmed the fragility: holding out `PAY0ZVBAI` eliminates the DNA-damaging vs Topo-I comparison entirely. The all-mechanism context-fixed residual range was nominal under within-context permutation (observed 0.1874; p=0.0343), but this did not survive the payload-level identifiability requirement because no exact context has replicated payload support in both mechanism arms.

**Interpretation**

The current data can show payload-specific differences within a few exact biological contexts, but they cannot distinguish a mechanism-level effect from payload identity. The observed DNA-vs-Topo pattern is scientifically interesting only as a descriptive payload/context pattern, not as evidence for mechanism-aware response prediction.

**Decision**

**PIVOT-WEAK**

There is an interesting descriptive pattern, but the current dataset cannot separate mechanism from payload identity strongly enough for a defensible mechanism-level claim.

**Next cheapest discriminating experiment**

Stop modeling. To resolve the pivot, obtain or curate crossed data where the same antigen x cell-line contexts are measured under at least two independent payloads per mechanism, across several contexts per mechanism pair, with harmonized endpoints and units.
