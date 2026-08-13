# Chapter 2. Finding real synthetic-lethal partners: a DepMap screen

*Chapter 1 ended with a warning: expression data alone cannot tell a lethal partnership from a regulatory one. This chapter builds the model that can, and watches it rediscover textbook biology from scratch.*

## Where Chapter 1 left us

The MMP14–PTEN story failed for an honest reason. Co-silencing avoidance, two genes almost never switched off together, is easy to find in expression data and impossible to interpret on its own, because regulation, lineage, and cell state all forge the same pattern. To find real synthetic lethality you have to watch what happens when a gene is actually removed. That data exists, genome-wide, in DepMap.

Before building the screen, it is worth understanding *which* synthetic lethals it will find most cleanly, because one class stands apart, and two of our three test cases belong to it.

## Why paralogs are the cleanest synthetic lethals

Over evolutionary time, genes are duplicated. The copies, called **paralogs**, often keep overlapping functions, and that redundancy is the point: it is a built-in backup. A cell can lose one paralog and barely notice, because its twin covers the same job. This is precisely why so many genes look **dispensable** in a single-gene knockout screen; their importance is masked by a sibling quietly doing the same work.

But redundancy is not the same as unimportance. If the *shared* function is essential, then removing **both** paralogs pulls away the backup and the job goes undone, and the cell dies. That is synthetic lethality in its most mechanically transparent form: each gene dispensable alone, lethal together, for the simple reason that they were doing the same essential thing.

This makes paralogs the **most predictable and most learnable** class of synthetic lethal pair. You can often guess a candidate from sequence similarity and shared expression alone, and paralog pairs are heavily enriched among all experimentally validated synthetic lethals. They also matter enormously in cancer. Tumors frequently inactivate one member of a paralog family, by mutation or deletion, and in doing so become **addicted to the surviving twin**. That twin is now an Achilles heel present in the cancer and absent in normal tissue, which is exactly the therapeutic window synthetic lethality promises. The SWI/SNF chromatin-remodeling complex is the textbook example: tumors that lose **SMARCA4** lean on its paralog **SMARCA2**, and tumors that lose **ARID1A** lean on **ARID1B**.

That is why paralogs are the ideal yardstick for a discovery method: the right answer is known, and the biology is clean. A screen that cannot surface SMARCA2 in SMARCA4-mutant cells is not to be trusted with anything harder. Non-paralog, pathway-level synthetic lethals, such as PTEN and PI3Kβ, are subtler and more valuable to discover, so we include one to show the method reaches past the easy case as well.

## The idea: a natural experiment

DepMap knocks out every gene, one at a time, in about 1,100 cancer cell lines and measures how much each knockout hurts. It never removes two genes at once, but it doesn't have to. Synthetic lethality makes a sharp prediction:

> If genes A and B are synthetic-lethal, then knocking out A should hurt **specifically in the cells that have already lost B.**

Those B-deficient lines are already sitting in the panel. So for any tumor suppressor B, we can ask, gene by gene across the genome: is the dependency on A stronger in B-mutant cells than in B-wild-type cells? The genes where the answer is a clear yes are the candidate partners.

## The method, and the two things that make it trustworthy

For each gene we compare its knockout effect in B-mutant versus B-wild-type lines. Two choices separate a believable screen from a misleading one.

- **Define "deficient" by mutation, not expression.** A damaging mutation is a clean genetic loss. Low expression is not; it drifts with lineage and regulation, exactly the confounds that sank Chapter 1.
- **Correct for lineage.** PTEN and ARID1A losses cluster in particular cancer types, so a raw comparison can mistake "this is a glioma gene" for "this is a PTEN partner." We regress lineage out of both the dependency profiles and the mutation status before testing, so what remains is the within-lineage effect.

## Does it work? Point it at genes with known answers

The honest test of a discovery method is whether it rediscovers what we already know. We ran it on three tumor suppressors, each with a textbook synthetic-lethal partner, and asked where that partner landed.

![Lineage-corrected SL screens recover known partners](../sl_v2_panels.png)

*Each panel is one tumor suppressor. Green points are candidate partners (more required when the gene is lost); the known partner is marked in red.*

- **SMARCA4 to SMARCA2**: rank **#1** genome-wide. The classic SWI/SNF paralog buffering relationship, top of the list.
- **ARID1A to ARID1B**: rank **#1**, and **WRN**, an independently established ARID1A synthetic lethality, appears right behind it, alongside the SWI/SNF paralog ARID2.
- **PTEN to PIK3CB (p110β)**: rank **#3**. PTEN-null cancers depending on PI3Kβ is well-documented clinical biology, and the screen found it blind.

Three genes, three known partners at or near the very top: two paralogs and one pathway dependency. The engine works.

## The rigor dividend

There is a quiet lesson buried in the PTEN result. An earlier, looser version of this screen, with PTEN defined by *low expression* and no lineage correction, had produced an appealing shortlist: redox genes (GLRX3, TXNRD1) and a cell-cycle kinase (CDK2), a tidy "PTEN-null cells are under oxidative and proliferative stress" story.

None of them survived. Once PTEN loss was defined genetically and lineage was removed, only **PIK3CB** remained. The redox story had been an artifact of how "PTEN-deficient" was defined and which lineages happened to carry low expression. The same discipline that dissolved MMP14–PTEN in Chapter 1 dissolved these too, and left the one real hit standing.

This is the whole value of the method. Not that it produces a list, but that a well-built version produces a list you can believe.

## What this is, and what it still isn't

We now have a **validated, generalizable synthetic-lethal discovery engine**. Give it any tumor suppressor and it returns a ranked, lineage-corrected, statistically controlled list of candidate partners, anchored by the reassurance that it recovers known biology at rank one.

It is also **Tier 1** of the larger program. The per-gene statistic it produces is exactly the feature a supervised classifier (trained on curated synthetic-lethal databases) or a graph neural network would consume to predict pairs it has never seen. This screen is both an answer and a feature generator.

The limits are worth stating plainly, and Chapter 1 named the deepest one. These are 2D survival screens. They cannot see synthetic lethality that only emerges in three-dimensional matrix, under metabolic stress, or in the tumor microenvironment, the very regimes where a gene like MMP14 does its work. Mutation-defined loss also misses deletion-driven and epigenetic silencing. And a screen, however clean, still only nominates; the endpoint is a combinatorial CRISPR experiment on the top candidates. But that is now a short, well-posed list, recovered by a method that has earned its credibility on genes whose answers we already knew.

---

*How it was done, briefly: DepMap 26Q1 CRISPR gene effect across 1,116 cell lines; deficiency defined by damaging mutation (24Q4); for each gene, a lineage-corrected partial correlation between knockout effect and the target's mutation status (Frisch–Waugh–Lovell residualization on OncotreeLineage), Benjamini–Hochberg FDR. Ranked tables: sl_PTEN_v2.csv, sl_SMARCA4_v2.csv, sl_ARID1A_v2.csv.*
