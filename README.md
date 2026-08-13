# Synthetic lethality from cancer dependency data

**Baja Bio · Bioinformatics**

A toolkit for discovering and prioritizing cancer synthetic-lethal drug targets from public
genome-scale data (DepMap CRISPR dependencies, copy number, expression, mutation; Ensembl paralogs;
GTEx; protein interactions). It implements two complementary target-discovery models and the
driver-versus-passenger discipline that keeps their predictions honest.

The central idea: cancers delete DNA in blocks, so a deleted gene often drags an innocent neighbor
(a *passenger*) down with it. When that passenger has an essential **paralog**, the surviving copy
becomes a selective, druggable vulnerability. The same logic, applied to two lost genes at once,
finds **higher-order** (third-gene) targets.

---

## What is here

Two models, one principle.

1. **Paralog synthetic-lethality predictor** (`src/paralog_model/`)
   Predicts, for any paralog pair, whether one becomes essential when its partner is lost. Trained on
   DepMap-derived labels with generalizable features (sequence identity, family size, co-dependency,
   co-expression, protein interaction, normal-tissue expression). Recovers validated targets
   (VPS4A, SMARCA4, MAGOH) blind. Includes the honest finding that a DepMap-derived label is entangled
   with single-gene essentiality, and that generalizable pair features win on independent validation.

2. **Higher-order (third-gene) model** (`src/higher_order/`)
   Screens two-gene-loss backgrounds for a third gene that becomes selectively essential, corrected for
   tissue of origin, and classifies each hit as a genuine three-way interaction or a passenger-driven
   one. Produces a ranked target catalogue and per-cancer-type applications.

3. **Driver-versus-passenger analysis** (`src/dependency_screens/driver_vs_passenger.py`)
   Decomposes any two-loss dependency into its single-loss parts, separating a true interaction from a
   co-deleted passenger.

4. **Collateral lethality and dosage-sensitive targets** (`src/collateral_lethality/`)
   The ENO1/ENO2 homozygous-deletion case, its limitation for heterozygous events, and the
   expression-based search that finds dosage-sensitive targets (e.g., ATP1A1 for 1p/19q loss).

---

## Repository layout

```
src/
  data_prep/            Extract and stream DepMap matrices; call loss-of-function
  dependency_screens/   Lineage-corrected differential dependency; third-gene screen; driver/passenger
  paralog_model/        Build labels (v1 expression, v2 copy-number+complex+GTEx), train, validate
  higher_order/         Systematic third-gene model; pancreatic screens; cancer applications
  collateral_lethality/ ENO1/ENO2; locate 1p/19q partners; expression-based target search; ATP1A1
  applications/         Per-cancer profiles; cross-cancer paralog-target map; conservation fetch
results/
  models/               Trained model (paralog_sl_model.pkl) + ranked prediction tables
  tables/               Per-screen output tables (CSV)
  figures/              All figures (PNG)
docs/
  DATA_SOURCES.md       Every external dataset, with download URLs / figshare IDs
  chapters/             Long-form write-ups of each result (Markdown)
```

## Models and outputs (`results/`)

| File | What it is |
|---|---|
| `models/paralog_sl_model.pkl` | Trained paralog-SL classifier (gradient-boosted trees) + feature list |
| `models/paralog_SL_predictions.csv` | Every human paralog pair, scored |
| `models/third_gene_targets.csv` | Ranked third-gene targets by two-loss background, with driver/passenger call |
| `tables/pancreatic_paralog_targets.csv` | Pancreatic passenger-paralog target portfolio |
| `tables/oligo_1p19q_expression_targets.csv` | Dosage-sensitive targets for 1p/19q loss |
| `tables/sl_{PTEN,SMARCA4,ARID1A}_v2.csv` | Lineage-corrected dependency screens |

## Pipeline

```mermaid
flowchart TD
  A[CRISPR dependencies<br/>DepMap] --> D[Call loss of function<br/>mutation, deletion, silencing]
  B[Copy number, expression, mutation<br/>DepMap] --> D
  C[Paralog pairs + identity<br/>Ensembl] --> E[Map paralogs and<br/>co-deleted neighbors]
  D --> E
  E --> F[Differential dependency<br/>tissue-corrected]
  F --> G[Driver vs passenger<br/>single-loss test]
  G --> H[Ranked targets +<br/>cross-cancer map]
```

## Install

```bash
python -m venv venv && source venv/bin/activate
pip install -r requirements.txt
```

## Getting the data

The scripts operate on public datasets that are **not redistributed here** (DepMap has its own terms of
use). See [`docs/DATA_SOURCES.md`](docs/DATA_SOURCES.md) for every source with exact download URLs and
figshare file IDs. Place the downloaded files in a working directory and adjust the input paths at the
top of each script (they currently point to local filenames such as `gene_effect_full.csv`).

## Running

A typical flow:

1. `src/data_prep/` — download and subset the DepMap matrices; call loss of function.
2. `src/dependency_screens/third_gene_screen.py` — screen a two-gene-loss background.
3. `src/dependency_screens/driver_vs_passenger.py` — decompose hits into driver vs passenger.
4. `src/higher_order/systematic_third_gene_model.py` — the full third-gene target catalogue.
5. `src/paralog_model/01..05` — build labels, train, validate the paralog predictor.
6. `src/applications/cross_cancer_paralog_targets.py` — the cross-cancer target map.
7. `src/collateral_lethality/` — the oligodendroglioma / dosage-sensitive analyses.

## Honest caveats

This toolkit generates **hypotheses**, not validated targets.

- Dependencies come from **cell lines**, not patient tumors, in 2D culture.
- "Loss" is often approximated from **expression or copy number**, a proxy for the true genetic event.
- A predictive label built from single-gene dependency is **entangled with essentiality**; independent
  labels (combinatorial double-knockout screens) are the correct upgrade. See
  `docs/chapters/chapter-06-paralog-model.md`.
- Higher-order "genuine three-way" calls are a **first-pass** classification; most two-loss dependencies
  are single-driver plus a co-deleted passenger.
- On-target toxicity matters: a nominated target still needs a modality and a therapeutic window
  (see the ATP1A1 / antisense-oligonucleotide discussion in the collateral-lethality write-up).

## License

MIT. See [`LICENSE`](LICENSE). External datasets retain their own licenses and terms of use.
