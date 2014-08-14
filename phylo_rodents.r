# code to build a phylogenetic tree for the 21 species at Portal
# use the tree from Bininda-Emonds et al. 2007. The delayed rise of present-day mammals. 
#                   Nature. 446: 507-512

library(ape)
library(geiger)
library(ggplot2)
library(picante)
library(PhyloOrchard)
#if PhyloOrchard not already installed, for mac use the terminal: svn checkout svn://scm.r-forge.r-project.org/svnroot/phyloorchard/
#                                       for windows use: install.packages("PhyloOrchard", repos="http://R-Forge.R-project.org")

wd = "/Users/sarah/Documents/GitHub/portal-rodent-dispersal"
setwd(wd)
load(".Rdata") #load saved Rdata

# read in species trait data from earlier analysis/research/etc.
traits=read.csv("traits.csv", row.names=1) #comes from mall in rodent_wrapper

# get phylo data from the published mammal tree
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

# as an example, distance between taxa
trx["Dipodomys_merriami", "Chaetodipus_baileyi"]

# Lets use a subset of the traits we are interested in and order by the tree tips
#TODO: Need to translate the traits into z-scores so they are on the same scale
traitDF = traits[,c(2,4,5,6,7,11,12,13,16,18,20)]
traitDF <- traitDF[tree.p$tip.label, ] 

# Standardize the matrix to correct for different units by subtracting means and dividing by sd
zscoret = apply(traitDF, 2, function(x) {
  y = (x - mean(x))/sd(x)
  return(y)
})
rownames(zscoret) <- rownames(traitDF)
zscoret = as.data.frame(zscoret)

# make a pairs plot to look at correlation structure of standardized data
pairs(zscoret[,c(1,4,9,11)], pch = 19)

# the correlation structure expected if traits evolve by Brownian motion 
# and fit a generalized least squares model assuming this correlation structure.
bmRodents <- corBrownian(phy=tree.p) 

bm.gls <- gls(meanabun ~ propyrs, correlation = bmRodents, data = zscoret) 
summary(bm.gls)

bm.gls <- gls(S ~ Psi, correlation = bmRodents, data = zscoret) 
summary(bm.gls)

bm.gls <- gls(modal_distance ~ propyrs, correlation = bmRodents, data = zscoret) 
summary(bm.gls)



# Again we will first build the correlation structure, this time assuming if  
# traits evolve as expected under the Ornstein-Uhlenbeck process with a variance-restraining 
# parameter, alpha. Ape automatically estimates the best fitting value of alpha for your data.
ouRodents <- corMartins(1,phy=tree.p) 
ou.gls<-gls(S ~ Psi, correlation=ouRodents, data=zscoret) 
summary(ou.gls)

