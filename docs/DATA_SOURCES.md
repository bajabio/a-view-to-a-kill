# Data sources

None of these datasets are redistributed in this repository. Download them from the sources below and
place them where the scripts expect them (adjust the paths at the top of each script). DepMap, GTEx,
Ensembl, BioPlex, and SynLethDB each retain their own license and terms of use.

## DepMap (Broad Institute)

Portal: https://depmap.org/portal/ · Files are hosted on figshare.

| Release | File | figshare file ID | Used for |
|---|---|---|---|
| 26Q1 (Chronos) | `gene_effect.csv` (CRISPR dependency, ~413 MB) | `67214582` | Dependency scores (the "essentiality") |
| 24Q4 Public | `Model.csv` | `51065297` | Cell-line lineage (OncotreeLineage) |
| 24Q4 Public | `OmicsExpressionProteinCodingGenesTPMLogp1.csv` (~507 MB) | `51065489` | Gene expression |
| 24Q4 Public | `OmicsSomaticMutationsMatrixDamaging.csv` (~148 MB) | `51065747` | Damaging-mutation calls |
| 24Q4 Public | `OmicsAbsoluteCNGene.csv` (~239 MB) | `51065303` | Copy number (deletions) |

Download pattern: `https://ndownloader.figshare.com/files/<FILE_ID>` (use `curl -sSL`).
The large matrices are rows = ModelID, columns = `SYMBOL (EntrezID)`. Use
`src/data_prep/extract_*` to stream-subset the columns you need.

## Ensembl paralogs (BioMart)

Human protein-coding paralog pairs with pairwise sequence identity. The main host was intermittently
returning 404; the **useast mirror** worked:

```
https://useast.ensembl.org/biomart/martservice?query=<url-encoded XML>
```

BioMart XML query (dataset `hsapiens_gene_ensembl`): filter `biotype = protein_coding` and
`with_hsapiens_paralog`; attributes `external_gene_name`, `hsapiens_paralog_associated_gene_name`,
`hsapiens_paralog_perc_id`, `hsapiens_paralog_perc_id_r1`. Saved locally as `paralogs.tsv`.

## GTEx v8 (normal-tissue expression)

Median TPM per gene across tissues:
```
https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz
```

## BioPlex (protein interactions)

Used as the co-complex / physical-interaction feature (substituting for CORUM, whose host was down):
```
https://bioplex.hms.harvard.edu/data/BioPlex_293T_Network_10K_Dec_2019.tsv
```

## SynLethDB (curated synthetic-lethal pairs)

Obtained via the KG4SL release (Wang et al., Bioinformatics 2021):
```
https://raw.githubusercontent.com/JieZheng-ShanghaiTech/KG4SL/main/data/sl_data
https://raw.githubusercontent.com/JieZheng-ShanghaiTech/KG4SL/main/data/dbid2name.csv
```
Used for supervised-label experiments and comparison (see chapter 3 write-up).

## mygene.info (annotation)

- Conservation proxy: PANTHER ortholog breadth via `pantherdb.ortholog.taxid`.
- Chromosome arm / cytoband: `map_location` field.
Endpoint: `https://mygene.info/v3/query` (POST, `scopes=symbol`, `fields=...`, `species=human`).

## Note on the "combinatorial dual-knockout screens"

The truly independent labels for paralog synthetic lethality (Thompson, Dede, CHyMErA, Parrish, and the
De Kegel & Ryan compilation, Cell Systems 2021) come from journal supplementary tables that are behind
publisher access controls and were not machine-downloadable in this environment. They are the correct
label source for a production model; see the chapter 6 write-up.
