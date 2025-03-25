# MMM
Code to run BiSSE on Carnivore Mammals to determine if past MMM champions are spread through a phylogeny of carnivores randomly, or grouped into one or a few lineages.

DISCLAIMER: BiSSEE is absolutely not designed to test this dumb theory, instead if you have a real binary trait that exists in some but not all species, BiSSE will test if the distribution is non-random and if lineages with the trait speciate more, and go extinct less (or the inverse).

MMM_biased.R is the R code to run BiSSEE on the Carnivora data. Lines 1-59 get you to the ANOVA p-value, lines after are for running Speciation and Extinction estimates, which are not appropiate for what we're trying to do for MMM but can be run on real life history data.

mammal_phylogeny.phy is a phylogeny of 244 mammals from Heuer et al 2024: https://doi.org/10.1038/s41598-024-74747-0

mammal_data_reduced.csv has a list of all species and the column morph can be modified to introduce a binary trait as 1-0, in this case, has this species been a previous MMM champion.
