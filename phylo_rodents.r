# code to build a phylogenetic tree for the 21 species at Portal
# use the tree from Bininda-Emonds et al. 2007. The delayed rise of present-day mamals. 
#                   Nature. 446: 507-512

library(picante)
library(PhyloOrchard)
#if PhyloOrchard not already installed, for mac use the terminal: svn checkout svn://scm.r-forge.r-project.org/svnroot/phyloorchard/

data(BinindaEmondsEtAl2007)
plot(BinindaEmondsEtAl2007[[1]], show.node.label = TRUE)

tree = BinindaEmondsEtAl2007[[1]]

# see prune.sample, normally this takes in a siteXspp matrix, have to fool it here
d <- matrix(nrow = 1, ncol = 21)
colnames(d) <- c("Dipodomys_merriami", "Dipodomys_ordii", "Dipodomys_spectabilis",
                 "Chaetodipus_penicillatus", "Chaetodipus_baileyi", "Perognathus_flavus",
                 "Peromyscus_maniculatus", "Peromyscus_eremicus", "Peromyscus_leucopus",
                 "Reithrodontomys_megalotis", "Reithrodontomys_montanus", "Reithrodontomys_fulvescens",
                 "Baiomys_taylori", "Chaetodipus_hispidus", "Chaetodipus_intermedius",
                 "Sigmodon_hispidus", "Sigmodon_fulviventer", "Sigmodon_ochrognathus",
                 "Neotoma_albigula", "Onychomys_torridus", "Onychomys_leucogaster")

# prune the tree
tree.p <- prune.sample(samp=d, phylo=tree)
plot(tree.p)

# get distance matrix for all species pairs
trx <- cophenetic(tree.p)

# distance between taxa
trx["Dipodomys_merriami", "Chaetodipus_baileyi"]

