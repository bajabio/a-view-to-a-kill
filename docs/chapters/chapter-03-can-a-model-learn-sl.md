# Chapter 3. Can a machine learn to predict synthetic lethality?

*A model trained on 35,000 known synthetic-lethal pairs looks impressive, until you test it the way that matters. Then it turns out to be recognizing famous genes, not lethal partnerships.*

## The promise

Chapter 2 built a discovery engine, but a narrow one. Point it at a gene and it finds that gene's partners; every new target starts from scratch. The tempting next step is to build something general. We already know thousands of synthetic-lethal pairs from decades of experiments. Why not train a computer to study all of them at once, learn what a lethal pair "looks like," and then predict brand-new pairs we have never tested?

There is plenty to learn from. Public databases hold tens of thousands of experimentally reported synthetic-lethal pairs. We took about 35,000 of them as our positive examples, kept the ones where both genes were measured in DepMap, and added an equal number of random gene pairs to serve as negatives. Then we described each pair with two kinds of information, both drawn from the earlier chapters:

- **The pair relationship**: how similarly the two genes behave when knocked out across the cell-line panel, the same co-dependency signal that powered Chapter 2. Genes that work in the same complex or pathway tend to track together.
- **Each gene on its own**: simple summaries of how essential and how variable each single gene is across cell lines.

We fed all of this to a standard machine-learning model and asked it to tell real synthetic-lethal pairs from random ones. Throughout, think of the model's score on a simple scale: **1.0 is a flawless predictor and 0.5 is a coin flip.**

## It looks like it works

Tested in the usual way, holding back a random slice of the pairs and asking the model to score them, it reached **0.83**. That is a respectable number, the kind that fills tables in synthetic-lethality prediction papers. If we stopped here we would report a success. We should not stop here.

## The trap: it recognizes genes, not partnerships

Synthetic lethality is a property of a *pair*, but the database is built from *genes*, and the same gene sits in hundreds of different pairs. PTEN alone appears in a great many. When we hold back a random slice of pairs, PTEN is in the training examples a hundred times and then again in the test. The model does not have to understand anything about the PTEN partnership. It only has to memorize that "pairs containing PTEN are usually labelled lethal." It is passing the test by recognizing the students, not by knowing the material.

The honest test is to hold back whole **genes**, not pairs, and then score only the pairs in which *both* genes are ones the model has never seen. That is the situation that actually matters, because the entire purpose of a predictor is to say something about genes you have not already screened.

Tested that way, the score falls from 0.83 to **0.71**. That drop is the price of the memorization the first test allowed, and it is the single mistake that quietly inflates much of the published literature.

![A synthetic-lethality predictor looks strong until tested honestly](../ch3_auc.png)

*The same model under four tests. Letting it see familiar genes flatters it; forcing it onto genuinely new genes is more honest; and stripping it down to one kind of information at a time (right) shows where its remaining skill actually comes from.*

## What is it actually learning?

A 0.71 is not nothing, so the fair question is what carries it. We rebuilt the model twice, each time giving it only one kind of information, and tested both on genes it had never seen.

- **Using only "how essential is each gene on its own": 0.71.** Almost the whole model.
- **Using only the pair relationship: 0.51.** A coin flip.

This is the uncomfortable result. On genuinely new genes, the information that describes an actual relationship between the two partners adds essentially nothing. What the model has really learned is which *individual* genes are essential and heavily studied, because those are exactly the genes that fill synthetic-lethality databases. It is a fame detector wearing the costume of a partnership predictor.

## The smoking gun

If more proof were needed, the model hands it over. Ask it to score the pair from Chapter 1:

> **MMP14 and PTEN: predicted chance of synthetic lethality, 98%.**

The model is almost certain about the very pair that two chapters of evidence, the expression pattern, the published prostate and pancreatic cancer biology, and the stratified CRISPR test, showed is *not* synthetic-lethal. And it is certain for the worst possible reason. PTEN is a famous, essential, heavily-annotated tumor suppressor, and the model has learned that famous essential genes make lethal pairs. The seductive false positive from Chapter 1 walks right back in the front door, now wearing a confidence score.

## Why it happens, and what to do instead

This is not a flaw you can tune away. It is built into the ingredients.

- **The examples are lopsided.** Databases are dominated by well-studied genes, so "is this a famous essential gene" becomes an easy shortcut to the right answer.
- **The negatives are not truly negative.** A random pair labelled "not lethal" has usually just never been tested, so the model partly learns what people have looked at, not what is real.
- **Context is thrown away.** A pair can be lethal in one tissue or genotype and not another, but the database pools everything together, blurring the very thing we want to predict.

None of this makes the approach useless. It means three practical rules. Always test on genes the model has never seen, or you are grading it on memory. Trust it only where the pair relationship genuinely carries information, above all for **paralogs**, the duplicated-gene pairs from Chapter 2, where family resemblance gives the model something real and specific to hold onto. And treat its output as a **ranked list of suggestions**, never as evidence, with the direct dependency test from Chapter 2 as the check that a suggestion actually behaves like a real vulnerability.

## The through-line

Three chapters, three increasingly powerful methods, and three different ways to fool yourself. A correlation in expression confused ordinary regulation with lethality. A CRISPR dependency screen had to be corrected for tissue of origin and separated from an expression artifact. And a machine-learning model, the most powerful tool of the three, quietly learned which genes are famous and handed back the original false positive with 98% confidence.

The tool was never the point. What separated a real finding from an artifact, in every chapter, was the same handful of habits: define what you are measuring honestly, test in a way the method cannot game, take it apart until you know which piece is doing the work, and confirm against mechanism. The computer narrows the search. It does not end it. The last word still belongs to the experiment, and the work in between is making sure the short list you hand the person at the bench is short for the right reasons.

---

*How it was done, briefly: positives were about 35,000 human synthetic-lethal pairs from the SynLethDB database (via the KG4SL release) whose two genes were both present in DepMap; negatives were an equal number of random gene pairs from the same set of genes. Each pair was described by its co-dependency (how correlated the two genes' knockout effects are across 1,208 cell lines in DepMap 26Q1) and by simple single-gene essentiality summaries. A gradient-boosted-tree model was scored two ways: holding out random pairs (which lets genes reappear) versus holding out whole genes (so test pairs are entirely unseen), plus one-feature-at-a-time tests under the stricter, gene-held-out scheme.*
