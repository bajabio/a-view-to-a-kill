# Chapter 6. A model that prioritizes paralog targets, and the label that makes or breaks it

*Can a machine, shown every paralog pair we know, learn to point at the next drug target? It can, but only once you stop letting it grade itself on a shadow of the answer.*

## The promise, and its limit

The screens of the earlier chapters are powerful but narrow. Point one at a gene loss and it ranks the whole genome for a lethal partner, but only when enough cell lines carry that loss. Rare losses, or combinations never sampled, are invisible to it. The obvious next step is a model that generalizes: train it on every paralog pair we already understand, teach it what a lethal pair looks like, and let it score pairs it has never seen. This chapter builds that model for the class where it is most likely to work, the paralogs, the duplicated genes that buffer each other and that gave us VPS4A in Chapter 5.

## Building it

The ingredients are all in hand. Every human paralog pair with its sequence similarity comes from Ensembl. Whether one paralog becomes essential when its partner is lost comes from the CRISPR dependency map. And each pair can be described by properties that do not require running the experiment: how similar the two genes are, how large their gene family is, whether they share a protein complex, whether they rise and fall together in expression, how broadly the target is expressed in normal tissue. We trained a standard model on all of it, scoring pairs on a scale where 1.0 is a perfect predictor and 0.5 is a coin flip.

It works. With the three textbook cases held out of training entirely, the model ranks each in the top fraction of a percent of all paralog pairs.

![The model works, but leans on the wrong thing](../paralog_model.png)

*Left: held out entirely, the model places VPS4A, SMARCA4, and MAGOH among the top 0.2% of all paralog pairs. Right: but stripped to one kind of information at a time, the target's own selectivity does almost the whole job, while sequence identity barely helps.*

## The catch, again

The score looked excellent, an AUC of 0.98. But the honest number for a rare event is not that; it is the average precision, which was only 0.16. And the deeper problem was the same one from Chapter 3. When we gave the model one kind of information at a time, the target's own essentiality profile, whether it is a gene that is needed in some cell lines and not others, reproduced nearly the entire model. Sequence identity, the property that actually describes the *pair*, added almost nothing. The model was less a synthetic-lethality predictor than a detector of selective dependencies that happen to have a paralog.

We tried to fix this the way most people would, with better and more data. We swapped the crude expression-based definition of gene loss for real copy-number deletions. We added protein-complex membership and normal-tissue expression from GTEx. The score inched up. The behavior did not change at all. Essentiality still did the whole job; sequence identity still sat idle.

## The diagnosis: it is the label, not the features

The reason is not the features. It is the label. Our definition of a synthetic-lethal pair was built from the CRISPR data itself: a pair counts as lethal if the target becomes more essential when its partner is lost. But that definition is, almost by construction, a statement about the target being a selective dependency. So a feature that measures selective dependency will always predict it, and no amount of sequence, complex, or expression information can dislodge something that is baked into the definition of the answer. A model can only be as honest as the label it is trained against, and ours was quietly grading the model on a shadow of itself.

## The resolution

The way to break the circle is to judge the model against synthetic lethals defined *independently* of the dependency map, pairs confirmed in the literature and in combinatorial double-knockout experiments. There are thirteen such gold-standard paralog SLs. We held every one of their genes out of training, then asked two models to rank them: one allowed only the target's essentiality, and one allowed only the generalizable pair features, sequence identity, shared complex, co-expression, normal-tissue expression, with no essentiality at all.

![Judged on independent SLs, the pair features win](../indep_validation.png)

*Each validated pair, ranked by the two models. Lower is better. The model that knows only the pair relationship, and nothing about the target's essentiality, recovers the independent synthetic lethals better overall, and rescues the ones the essentiality model gets badly wrong.*

The result reverses the earlier picture. Judged against real, independent synthetic lethals, the model that uses *only* the pair relationship, and is told nothing about whether the target is a selective dependency, does better, a median rank in the top 0.6 percent against the top 1.3 percent for the essentiality model. More tellingly, it rescues the pairs the essentiality model fails outright. ME2 and ME3, metabolic paralogs whose target is not a dramatic stand-alone dependency, sat near the very bottom of the essentiality ranking and near the top of the pair-feature one. The same held for the SMARCC1 and SMARCC2 chromatin paralogs.

So the pair features were never useless. They were being outvoted by a label that rewarded the wrong thing. Given an honest label, sequence identity and shared complex membership are exactly what find real synthetic lethals, including the ones a single-gene view misses.

## The lesson

There is a hierarchy of ways to fool yourself, and this chapter reaches the deepest one. Chapter 1 warned that a pattern can have the wrong cause. Chapter 3 warned that a model can be graded so that it cheats. This chapter warns that even a well-evaluated model can be quietly wrong if the *label*, the definition of the thing you are predicting, is entangled with a feature you already have. You can hold out genes correctly, ablate honestly, and still be measuring a shadow, because the error is upstream of all of it, in how you defined success.

The practical answer is simple to state and easy to skip: define the target variable from a source independent of your features. For paralog synthetic lethality that means labels from combinatorial double-knockout screens, not from the same single-gene dependency data that supplies the features. Do that, and the model is no longer a fame detector or a selectivity detector. It becomes what it was meant to be, a predictor of the relationship between two genes, and it will hand you a ranked list of paralog targets that is worth taking to the bench.

Machine learning can prioritize real targets. Whether it does comes down, in the end, not to the model or even the features, but to whether the label measures the thing you actually care about.

---

*How it was done, briefly: human paralog pairs with sequence identity from Ensembl BioMart; per-pair labels from DepMap 26Q1 CRISPR (target more essential in partner-lost lines, with loss defined by copy-number deletion or low expression); features included sequence identity, family size, co-dependency, co-expression, protein interaction from BioPlex, and GTEx normal-tissue expression. Models were gradient-boosted trees evaluated with genes held out of training. Independent validation used thirteen literature-confirmed paralog synthetic lethals, with all of their genes excluded from training.*
