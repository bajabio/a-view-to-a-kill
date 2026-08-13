# Chapter 4. Higher-order synthetic lethality: finding the third gene

*When a cancer has already lost two genes, is there a third whose knockdown finishes it off? A genome-wide screen finds real answers, and a trap that looks just like one.*

## Beyond pairs

Everything so far has been about pairs: gene A and gene B, one lethal partnership at a time. But real tumors are not so tidy. A single cancer genome routinely carries several inactivated genes at once, and those losses do not act in isolation. If losing one tumor suppressor can create a dependency, as Chapter 2 showed, then losing two together might reshape the cell more drastically still, opening a vulnerability that neither loss creates on its own. The question is whether a third gene, knocked down on top of a specific two-gene loss, becomes the thing the cell can no longer survive without.

This is higher-order synthetic lethality, and it is exactly the kind of thing a dependency map is built to test. The logic is the same as before, only the background is richer. Instead of asking which gene becomes essential when a single partner is lost, we ask which gene becomes essential in the cells that have lost a particular *pair*.

## The method

We first had to define loss honestly. A gene can be knocked out of action by a damaging mutation or by outright deletion, and the two leave different fingerprints, so we called a gene lost if it carried a damaging mutation or if its expression had collapsed into the bottom of the range, which catches deletions and silencing. For every pair of tumor suppressors of interest, we then found the cell lines that had lost both, and there were enough of them to work with: a hundred or so cell lines carried the classic ninth-chromosome co-deletion, and dozens carried each of the other combinations.

With a double-loss background defined, the screen is Chapter 2's engine pointed at a richer target. Across the whole genome, we ask whether each gene is more essential in the double-loss cells than in the rest, after correcting for tissue of origin so that a vulnerability of one cancer type is not mistaken for a consequence of the two losses.

## It works: a textbook vulnerability, recovered blind

The cleanest test is a background whose answer is already known. Tumors that delete the *CDKN2A* locus very often delete the neighboring gene *MTAP* at the same time, because the two sit side by side on the chromosome. *MTAP* loss is famous for creating a dependency on the arginine methyltransferase *PRMT5*, a vulnerability now being pursued in the clinic.

Run blind on the *MTAP* and *CDKN2A* double-loss background, the screen recovers that biology in full.

![Two-gene-loss to third-gene lethality screen](../thirdgene_result.png)

*Left: in cells that have lost MTAP and CDKN2A, the genes that become most selectively essential are the PRMT5 methylation axis, with WDR77 (PRMT5's partner protein) at the top, alongside PRMT5 itself, PELO, PRPF6, and METTL16. Middle: cells that have lost TP53 and RB1 become dependent on the E2F and cell-cycle machinery. Right: the caveat, explained below.*

The top of the list is not PRMT5 in isolation but the whole methylation module: **WDR77**, the partner protein that PRMT5 cannot work without, ranks first, with **PRMT5** just behind it, joined by **PELO**, **PRPF6**, and **METTL16**, all genes tied to the same axis. A method that reconstructs an entire drug-target pathway from nothing but a genetic background is a method worth trusting with an unknown one.

## A real discovery: losing both brakes

The unknown one is compelling. Cells that have lost both *TP53* and *RB1*, the two master brakes on the cell cycle, become dependent on the cell-cycle engine itself. The strongest hit is **E2F3**, the transcription factor that RB1 normally restrains, followed by **SKP2**, **CKS1B**, **CDK2**, the cyclins **CCNE1** and **CCNE2**, and the replication-licensing genes **CDT1** and **ORC6**. The interpretation is almost mechanical: with both restraints gone, the cell races through division and can no longer tolerate losing any part of the machinery that drives it. That is a coherent, druggable dependency created specifically by the pair of losses. A third background, the co-mutation of *STK11* and *KEAP1* in lung cancer, pointed more weakly toward energy and mitochondrial genes, a suggestive lead rather than a firm one.

## The trap: a passenger dressed as an interaction

Here is where the discipline of the earlier chapters earns its keep. It is tempting to read "PRMT5 is essential in *MTAP* and *CDKN2A* double-loss cells" as a three-way interaction, a vulnerability that needs both losses. It is not. When we break the cells into groups, PRMT5 is just as essential in cells that have lost *MTAP* alone as in the double-loss cells, and cells that have lost only *CDKN2A* are no more dependent on PRMT5 than normal cells are. The dependency is driven entirely by *MTAP*. *CDKN2A* is only along for the ride, deleted because it happens to sit next to *MTAP* on the chromosome, a passenger mistaken for a partner.

This matters because it is the default outcome, not the exception. Across the backgrounds we tested, most of the third-gene dependencies were carried by whichever of the two losses was the stronger driver, with the second loss contributing little. A genuine higher-order interaction, a gene that needs both losses and neither alone, does exist, but it is quantitatively modest: the clearest recurrent examples were DNA-repair and mitotic genes such as MMS22L, CDK2, SKP2, and KIF18A in *TP53*-loss backgrounds, and the extra lethality from the second hit was small. The E2F dependency in the *TP53* and *RB1* background is likewise led by the loss of *RB1*, the direct brake on E2F, with *TP53* loss deepening rather than creating it.

## What this adds

Two-gene-loss backgrounds are real, they are common in cancer, and they do create third-gene dependencies that a screen can find and that sometimes point to genuine drug targets. But a double-loss background is very often, underneath, a single-driver dependency with a passenger attached, and the interesting biological claim, that a vulnerability requires *both* hits, is the harder one to earn. Distinguishing the two takes exactly the move this chapter made at the end: not stopping at "essential in the double-loss cells," but asking whether either loss alone already explains it.

It is the same lesson the series keeps arriving at from new directions. A pattern that looks like a specific interaction, a co-silencing, a co-dependency, a context vulnerability, is worth only as much as the test that separates the real interaction from the thing that merely resembles one.

---

*How it was done, briefly: loss of function was called from DepMap 24Q4 data as a damaging mutation or expression below the fifteenth percentile per gene. For each two-gene-loss background, genes more essential in the double-loss cell lines were identified from DepMap 26Q1 CRISPR data by a lineage-corrected partial correlation between the double-loss status and each gene's knockout effect. Whether a dependency required both losses was assessed by comparing knockout effects across double-loss, single-loss, and wild-type groups.*
