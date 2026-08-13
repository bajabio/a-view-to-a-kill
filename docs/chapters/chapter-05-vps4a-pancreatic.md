# Chapter 5. A pancreatic-cancer target from a passenger deletion: VPS4A

*Point the higher-order screen at the real genetics of pancreatic cancer, and it finds a drug target. The target does not come from the tumor suppressor that was lost, but from the innocent gene deleted next to it.*

## Real genetics, not a toy background

The previous chapter built a screen for a lethal third gene given two lost ones, and tested it on backgrounds chosen for convenience. The honest test is a real cancer, and pancreatic ductal adenocarcinoma is almost designed for the purpose. Its genetics are stereotyped: an activating mutation in KRAS starts the disease, and then three tumor suppressors are lost, in a rough order, as it progresses. CDKN2A goes first, in around sixty percent of tumors, then TP53 and SMAD4, each in roughly half. A pancreatic cancer is, more often than not, a cell that has lost two or three of these at once. So the question of the last chapter becomes concrete: given the losses that actually define this disease, what third gene can no longer be spared?

## The partners that did not hold up

The literature already offers candidates for SMAD4-loss in particular. Targeted CRISPR screens have nominated RAB10, KLF5, AURKA, and the BET proteins as synthetic-lethal partners of SMAD4. It is worth asking whether these behave, across the full panel of cancer cell lines, as genotype-selective dependencies, the way a real synthetic lethal should.

Mostly, they do not. When we test whether each becomes more essential in SMAD4-lost cells across the dependency map, RAB10 sits near the middle of the genome, AURKA is broadly essential everywhere and no more so when SMAD4 is gone, and the BET genes if anything lean the other way. Only KLF5 shows a weak selective signal. This is not a claim that the original papers are wrong; those interactions were found and validated in particular cell lines and contexts. It is a reminder that a hit from one screen is not yet a general dependency, and that the gap between the two is exactly what a pan-cancer test is for.

## The target the screen found on its own

Run blind on cells that have lost both SMAD4 and TP53, the strongest and most statistically secure hit is a gene the literature above never mentions: **VPS4A**. It sits at the top of the list, comfortably clear of the multiple-testing threshold, a gene that becomes selectively essential exactly in this double-loss background.

The reason is a small masterpiece of cancer genome geography. VPS4A has a paralog, **VPS4B**, a near-twin that does the same job in the cell, the two together running the ESCRT machinery that pinches off membranes during cell division. VPS4B happens to lie on the long arm of chromosome 18, about twelve million bases from SMAD4. When a pancreatic cancer deletes the piece of 18q that carries SMAD4, it very often takes VPS4B along with it, and that co-deletion is common: roughly two thirds of pancreatic cancers lose VPS4B this way. A cell that has lost VPS4B still has VPS4A, and gets by. Take away VPS4A as well, and there is no backup left for an essential job, and the cell dies. It is the paralog synthetic lethality of Chapter 2, delivered by the passenger biology of Chapter 4.

![VPS4A is a synthetic-lethal target in SMAD4-deleted pancreatic cancer](../pdac_vps4a.png)

*Left: in cells that have lost SMAD4 and TP53, VPS4A is the top selective dependency. Right: VPS4A becomes essential precisely when its paralog VPS4B is lost, the 18q deletion that also removes SMAD4; pancreatic cell lines, dark, are enriched in the vulnerable range.*

The data tell the mechanism directly. Sort the cell lines by how much VPS4B they still express, and VPS4A dependence tracks it cleanly: cells with the lowest VPS4B are the ones that cannot survive losing VPS4A, a relationship that is about as statistically firm as these data produce. Pancreatic lines are enriched among the VPS4B-low, vulnerable cells, as the co-deletion frequency predicts. This is not a subtle context effect that needs both losses in a particular combination; it is a clean, single-driver dependency created by the deletion of one gene, and VPS4A inhibitors are being pursued in the clinic for exactly this population.

## The lesson: look at the neighborhood, not just the driver

There is a pattern here that is more useful than any single target. Twice now, the most robust third-gene vulnerability in a cancer-defining loss has come not from the biology of the lost tumor suppressor but from an innocent bystander deleted alongside it. SMAD4 loss on 18q drags down VPS4B and creates a need for VPS4A. CDKN2A loss on chromosome 9 drags down the neighboring gene MTAP and creates a need for the PRMT5 methylation axis, which the same screen surfaces in the CDKN2A-lost background. In both cases the tumor suppressor is the marker, and the passenger is the mechanism.

This is good news for finding targets. A dependency that rides on a copy-number deletion comes with its own companion diagnostic, the deletion itself, and it points at a specific, druggable protein with a specific reason to be essential. The move that finds it is not to ask what the lost tumor suppressor does, but to look at what sits next to it on the chromosome, and to ask whether anything the cell needs was quietly lost in the same stroke.

## The through-line

The screen did its job here because it was anchored on the real genetics of a real disease, and because the final step was the one the whole series keeps returning to: when a gene came up as essential in a genotype, we did not stop at the association, but asked what actually drove it, and found a specific paralog loss rather than a vague consequence of the two mutations. The association was a starting point. The mechanism, a deleted paralog leaving its twin holding an essential job alone, is the thing you can take to an experiment, and in this case someone already has.

---

*How it was done, briefly: pancreatic driver losses were taken from the literature; loss of function was called from DepMap 24Q4 data as a damaging mutation or expression below the fifteenth percentile. Genes selectively essential in the SMAD4 and TP53 double-loss background were identified from DepMap 26Q1 CRISPR data by a lineage-corrected partial correlation. The VPS4A and VPS4B relationship was confirmed by comparing VPS4A knockout effect against VPS4B expression across 1,116 cell lines. Reported SMAD4 partners were assessed for genotype-selective dependency by the same lineage-corrected test.*
