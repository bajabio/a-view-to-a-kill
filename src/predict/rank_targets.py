# Rank-order third-gene targets for a two-gene-loss background, tissue-aware.
#
# Baja Bio · Bioinformatics — synthetic-lethality target-discovery toolkit
# Primary author: Jeff Milton
# Public dependency/genomic data per docs/DATA_SOURCES.md; set input paths for your setup.
"""
Given two background genes (A, B) and an optional tissue, return a rank-ordered list of candidate
third-gene targets: genes that become selectively essential when both A and B are lost.

The score is a lineage-corrected differential dependency (more negative = more selectively essential in
the double-loss background). Every hit is decomposed into its single-loss parts so a genuine
higher-order interaction is separated from a co-deleted passenger.

Loss is called from a boolean loss-of-function matrix (mutation OR deletion OR silencing) when one is
supplied, and falls back to expression (bottom `loss_pct` percentile) otherwise. The matrix matters:
a gene like TP53 is lost mainly by point mutation, which low expression does not see, so an
expression-only call misses most TP53-null lines and can return nothing.

Example
-------
    from rank_targets import load_data, rank_targets
    data = load_data("gene_effect_full.csv", "uni_expr.csv", "Model.csv",
                     loss_matrix_csv="lof_matrix.csv")
    df = rank_targets("SMAD4", "TP53", tissue="Pancreas", **data)
    df.head(20)
"""
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


def load_data(gene_effect_csv, expression_csv, model_csv, loss_matrix_csv=None):
    """Load and align the DepMap inputs. Returns a dict to splat into rank_targets().

    loss_matrix_csv (optional): a boolean ModelID x gene table of loss-of-function calls
    (mutation OR deletion OR silencing). Preferred over expression for genes lost by mutation
    (e.g. TP53), which low expression does not capture. Build it with data_prep/call_loss_of_function.py.
    """
    G = pd.read_csv(gene_effect_csv, index_col=0)
    G.columns = [c.split(" (")[0] for c in G.columns]
    E = pd.read_csv(expression_csv, index_col=0)
    lin = pd.read_csv(model_csv).set_index("ModelID")["OncotreeLineage"]
    lines = G.index.intersection(E.index)
    L = None
    if loss_matrix_csv is not None:
        L = pd.read_csv(loss_matrix_csv, index_col=0).reindex(lines).fillna(False).astype(bool)
    return {"G": G.loc[lines], "E": E.loc[lines], "lin": lin.reindex(lines), "L": L}


def rank_targets(gene_a, gene_b, G, E, lin, L=None, tissue=None,
                 loss_pct=15, min_double=15, fdr_max=0.25):
    """Rank candidate third-gene targets for the (gene_a, gene_b) double-loss background.

    Parameters
    ----------
    gene_a, gene_b : str      background genes (must be present in the expression matrix)
    G, E, lin      : from load_data()
    tissue         : str or None   OncotreeLineage to spotlight (e.g. "Breast", "Pancreas")
    loss_pct       : int      a gene is 'lost' below this expression percentile
    min_double     : int      minimum double-loss cell lines to proceed
    fdr_max        : float    keep hits at or below this Benjamini-Hochberg FDR

    Returns a DataFrame ranked by selective dependency, with a driver/passenger interpretation and,
    if a tissue is given, the effect within that tissue.
    """
    lines = G.index
    genes = np.array(G.columns)
    Gv = G.values.astype(np.float32)
    cm = np.nanmean(Gv, 0); ii = np.where(np.isnan(Gv)); Gv[ii] = np.take(cm, ii[1])

    def _loss(g):
        # prefer a mutation/deletion/silencing call from the loss matrix; else expression low
        if L is not None and g in L.columns:
            return L[g].values
        if g in E.columns:
            return (E[g] < np.nanpercentile(E[g], loss_pct)).values
        raise ValueError(f"no loss call available for {g} (not in loss matrix or expression)")
    lowA = _loss(gene_a)
    lowB = _loss(gene_b)
    both = lowA & lowB
    Aonly = lowA & ~lowB
    Bonly = ~lowA & lowB
    n = int(both.sum())
    if n < min_double:
        raise ValueError(f"only {n} double-loss lines for {gene_a}+{gene_b}; need >= {min_double}")

    # lineage-corrected differential dependency (Frisch-Waugh-Lovell)
    D = pd.get_dummies(lin, dummy_na=True).astype(float); D.insert(0, "_int", 1.0)
    Dv = D.values
    beta, *_ = np.linalg.lstsq(Dv, Gv, rcond=None)
    Xr = Gv - Dv @ beta
    y = both.astype(float)
    yr = y - Dv @ np.linalg.lstsq(Dv, y, rcond=None)[0]
    r = (Xr.T @ yr) / np.sqrt((Xr ** 2).sum(0) * (yr ** 2).sum() + 1e-12)
    dfree = len(lines) - np.linalg.matrix_rank(Dv) - 1
    t = r * np.sqrt(dfree / np.clip(1 - r ** 2, 1e-9, None))
    p = 2 * stats.t.sf(np.abs(t), dfree)
    fdr = multipletests(p, method="fdr_bh")[1]

    eff_double = Gv[both].mean(0)
    eff_a = Gv[Aonly].mean(0) if Aonly.sum() >= 5 else np.full(len(genes), np.nan)
    eff_b = Gv[Bonly].mean(0) if Bonly.sum() >= 5 else np.full(len(genes), np.nan)
    worse_single = np.nanmin(np.vstack([eff_a, eff_b]), axis=0)
    synergy = eff_double - worse_single  # < 0 => double worse than either single

    tis_mask = (lin == tissue).values if tissue else None
    eff_tissue = np.full(len(genes), np.nan)
    if tissue is not None:
        m = tis_mask & both
        if m.sum() >= 3:
            eff_tissue = Gv[m].mean(0)

    def interp(i):
        if not (t[i] < -3 and eff_double[i] < -0.4):
            return "weak"
        if synergy[i] < -0.1 and t[i] < -3.5:
            return "genuine higher-order"
        drv = gene_a if (not np.isnan(eff_a[i]) and eff_a[i] <= np.nanmin([eff_b[i], 0])) else gene_b
        return f"driven by {drv}"

    out = pd.DataFrame({
        "target": genes, "t": t, "fdr": fdr,
        "eff_double": eff_double, "eff_A_only": eff_a, "eff_B_only": eff_b,
        "synergy": synergy, "eff_in_tissue": eff_tissue,
        "interpretation": [interp(i) for i in range(len(genes))],
    })
    out = out[out.target.isin([gene_a, gene_b]) == False]
    out = out[(out.t < 0) & (out.eff_double < -0.4) & (out.fdr <= fdr_max)]
    out = out.sort_values("t").reset_index(drop=True)
    out.attrs["n_double_loss"] = n
    out.attrs["background"] = f"{gene_a}+{gene_b}"
    out.attrs["tissue"] = tissue
    return out


def main(argv=None):
    import argparse, sys
    ap = argparse.ArgumentParser(
        prog="a_view_to_a_kill",
        description="Rank third-gene targets for a two-gene-loss background "
                    "(add two genes as the cancer background, get a rank-ordered kill list).")
    ap.add_argument("gene_a"); ap.add_argument("gene_b")
    ap.add_argument("--tissue", default=None)
    ap.add_argument("--gene-effect", default="gene_effect_full.csv")
    ap.add_argument("--expression", default="uni_expr.csv")
    ap.add_argument("--model", default="Model.csv")
    ap.add_argument("--loss-matrix", default=None,
                    help="boolean ModelID x gene loss-of-function table; preferred for mutation-lost "
                         "genes such as TP53 (falls back to expression when omitted)")
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)
    data = load_data(a.gene_effect, a.expression, a.model, loss_matrix_csv=a.loss_matrix)
    df = rank_targets(a.gene_a, a.gene_b, tissue=a.tissue, **data)
    print(f"# {df.attrs['background']}  tissue={df.attrs['tissue']}  "
          f"double-loss lines={df.attrs['n_double_loss']}", file=sys.stderr)
    (df.to_csv(a.out, index=False) if a.out else df.head(30).to_csv(sys.stdout, index=False))


if __name__ == "__main__":
    main()
