# for manipulating  rodent movement data in MARK
## analyze survival and dispseral probabilities for species and guilds

library(RMark) #can't install on mac? Update R version?

#mac access
wd = "/Users/sarah/Desktop/portal-rodent-dispersal/"
setwd(wd)

source("movement_fxns.r") #functions for stake_movement.r
source("stake_movement.r") #makes a mark data structure using species-level data from Portal Project


# bring in the inp file and convert it to RMark format (.txt needs to be .inp, change manually)
MSdata <- convert.inp("do_mark.inp",
                      group.df=data.frame(sex=c("male","female","unidsex")))

# Build up the model. Looking at sex effects on dispersal/survival
MS.process <- process.data(MSdata,model="Multistrata",begin.time=2000,
                           group=c("sex"))

MS.ddl <- make.design.data(MS.process)

# Add dynamic dummy variable age class fields to the design data for Psi and p
MS.ddl$Psi$hy=0
MS.ddl$Psi$hy[MS.ddl$Psi$sex==0&MS.ddl$Psi$stratum=="A"]=1
MS.ddl$Psi$ahy=0
MS.ddl$Psi$ahy[MS.ddl$Psi$sex>=1&MS.ddl$Psi$stratum=="A"]=1

MS.ddl$p$hy=0
MS.ddl$p$hy[MS.ddl$p$sex==1&MS.ddl$p$stratum=="A"]=1
MS.ddl$p$hy[MS.ddl$p$sex==1&MS.ddl$p$stratum=="B"]=1
MS.ddl$p$sy=0
MS.ddl$p$sy[MS.ddl$p$sex==2&MS.ddl$p$stratum=="A"]=1
MS.ddl$p$asy=0
MS.ddl$p$asy[MS.ddl$p$sex>=3&MS.ddl$p$stratum=="A"]=1
MS.ddl$p$ahy=0
MS.ddl$p$ahy[MS.ddl$p$sex>=2&MS.ddl$p$stratum=="B"]=1

##### Add dummy variables for operating on specific states or transitions
  # A = 1 (home), B = 2 (away)
  # moving to B from anywhere, is risky
  # moving to A from anywhere, is less risky (within home, "normal" movements)
MS.ddl$Psi$toB=0
MS.ddl$Psi$toB[MS.ddl$Psi$stratum=="1" & MS.ddl$Psi$tostratum=="2"]=1

MS.ddl$Psi$toA=0
MS.ddl$Psi$toA[MS.ddl$Psi$stratum=="2"]=1

MS.ddl$p$strA=0
MS.ddl$p$strA[MS.ddl$p$stratum=="A"]=1
MS.ddl$p$strB=0
MS.ddl$p$strB[MS.ddl$p$stratum=="B"]=1

## TODO: fix recapture probabilities for unsampled or omitted months!
    # periods = c(267, 277, 278, 283, 284, 300, 311, 313, 314, 318, 321, 323, 
              #   337, 339, 344, 351)

# Add a field for monthly reporting probabilities 
MS.ddl$p$rpt=0
MS.ddl$p$rpt[MS.ddl$p$prd] == 267 = 0
MS.ddl$p$rpt[MS.ddl$p$prd] == 277 = 0

# TODO: Lots of work on building up the models!
                          