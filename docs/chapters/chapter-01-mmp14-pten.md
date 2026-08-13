# Chapter 1. Anatomy of a synthetic-lethal candidate: the MMP14–PTEN case

*A candidate synthetic-lethal pair, spotted in expression data, that every line of evidence we can bring to bear points away from, and toward a familiar PI3K/AKT regulatory axis, while leaving the definitive double-knockout still undone.*

## A tempting shortcut

Every cancer biologist knows why synthetic lethality is worth chasing. Two genes, each dispensable on its own but lethal to lose together; a tumor that has already dropped one partner is left leaning on the other, and that survivor becomes a target you can hit without harming normal tissue. BRCA and PARP made the principle famous, and the hunt for the next pair is a serious enterprise.

The gold standard for finding such pairs is a CRISPR or RNAi screen, powerful but expensive. So there is a standing temptation to nominate candidates first from data you already have. Expression data is the obvious place to look, because synthetic lethality should leave a mark there. If losing both genes kills a cell, then among the cells that survive you should almost never find both genes switched off at once, even though each is silenced on its own in plenty of other samples. The pair avoids being lost together.

We went looking for exactly that pattern across a large collection of expression measurements spanning hundreds of cell lines and tissues.

## A clean signal

One pair stood out: **MMP14 and PTEN**.

Each gene is silenced on its own in a good number of samples. Yet the two are almost never silenced in the same sample, far less often than chance would allow. Where you would expect a handful of doubly-silenced cells, there was essentially one.

![Candidate pairs: co-silencing is depleted](../sl_candidates.png)

*Each point is a sample; both axes run from expressed to silenced. The upper-right corner, where both genes are off, is where doubly-silenced survivors would sit, and for MMP14–PTEN (top left) it is empty.*

On its face, a promising lead. PTEN is a bona fide tumor suppressor, the pair is biologically plausible, and the pattern held up when we stripped out the usual confounds (genes that simply mark different cell lineages, and the broad proliferation-versus-quiescence program that makes many genes look anti-correlated). It looked like a synthetic-lethal hypothesis worth taking seriously.

## The catch

Here is the problem, and it is a conceptual one, not a statistical artifact. **Two very different biologies produce the identical pattern.**

- *Synthetic lethality:* cells die if both genes are lost, so the doubly-silenced state is missing because those cells are dead.
- *Regulation:* one gene actively controls the other, so losing the first pushes the second up, and the doubly-silenced state is missing because it never forms.

An expression pattern of "never both off" cannot, by itself, tell these apart. You need orthogonal evidence.

## What was already known

The literature settled the first half of the question immediately. In **PTEN-null prostate cancer cells**, loss of PTEN drives *more* MT1-MMP (MMP14) to the cell surface, and that increase runs squarely through **PI3K/AKT** signaling. The same inverse coupling has been reported in pancreatic cancer. So the relationship we picked up is real, but it looks regulatory. Losing PTEN doesn't obviously make a cell that also lacks MMP14 die; it turns MMP14 *on*. And tellingly, MMP14 has never appeared among PTEN's genuine synthetic-lethal partners (those are chromatin and replication factors such as BET, p300/CBP, and WDHD1).

That pointed away from synthetic lethality. But there was a way to test the lethality claim more directly.

## The test DepMap can actually do

Here we have to be careful. DepMap knocks out one gene at a time; it never removes both at once, so it cannot directly show that losing MMP14 *and* PTEN together is lethal. What it can do is lean on a natural experiment. Synthetic lethality predicts that MMP14 should become essential **specifically in cells that have already lost PTEN**. Those PTEN-deficient lines are already in the panel, so a single-gene MMP14 knockout within them delivers the second hit for free.

So we compared MMP14 knockout across about 1,100 cell lines, split by how much PTEN each one expresses.

![MMP14 dependency versus PTEN status](../sl_stratified.png)

*MMP14 knockout effect against PTEN expression. If the pair were synthetic-lethal, the PTEN-low cells on the left should sink below the dependency line; losing MMP14 should start to hurt once PTEN is gone.*

None of them do. In the most PTEN-deficient lines, MMP14 knockout has essentially no effect (the strongest response is a negligible −0.08), and there is no trend linking lower PTEN to greater MMP14 dependence. PTEN itself behaves exactly as a tumor suppressor should; losing it tends to *help* growth, not hurt it.

This is not a controlled double knockout, and it cannot be a definitive one. CRISPR knockout is rarely complete, MT1-MMP has redundant relatives among the other MMPs, and much of its biology plays out on the extracellular matrix, which a survival screen in a dish barely registers. But it is the most direct genetic test available, and it points the same way as everything else.

## The takeaway

No one has removed both genes together, so no single experiment proves the pair viable or lethal outright, and it is worth being honest that the direct double-knockout has not been done. But every line of evidence we *can* bring to bear points the same way. The expression pattern, the published in vivo biology, and the CRISPR dependency data all agree that MMP14 and PTEN are genuinely connected, through PI3K/AKT, as regulator and target, and none of them support a lethal partnership.

The broader lesson is one worth keeping close: **anti-correlation in expression is not genetic dependency.** Co-silencing avoidance is easy to find and easy to over-read, because regulation, lineage identity, and cell state all mimic the synthetic-lethal pattern. Used honestly, as a fast, cheap way to nominate candidates for a real screen, the approach earns its place. Used as evidence on its own, it will hand you regulatory relationships dressed up as drug targets.

The signal was real. The synthetic-lethal reading was not supported. Knowing the difference, and knowing which experiment would settle it, is the whole job. Here, that experiment is a clean one: knock out MMP14 in an isogenic PTEN-wild-type / PTEN-null pair, or run a dual-guide combinatorial screen, and the question is answered directly.

---

*How it was done, briefly: expression was summarized per gene per cell line and each gene called on or off; pairs were flagged when the both-off state was rarer than chance, after filtering out lineage-marker and cell-state effects. The synthetic-lethal claim was then tested in DepMap (26Q1 CRISPR gene effect; 24Q4 expression) by asking whether MMP14 knockout becomes damaging in PTEN-low cell lines; across 1,116 lines, it does not.*
