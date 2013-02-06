# for manipulating  rodent movement data in MARK
## analyze survival and dispseral probabilities for species and guilds
# Use multistrata models that give S(survival), p(capture probability) and Psi(transition probability)
# Will compare data among species and guilds. 
# Major covariates include sex, guild, and average body mass

library(RMark) 

#---------------------------------------------------------------------------------
#          bring in the data and source files
#---------------------------------------------------------------------------------
#set working directory and import source code
wd = "C://Users//sarah//Documents//GitHub//portal-rodent-dispersal"
setwd(wd)

source("stake_movement.r") #makes a mark data structure using species-level data from Portal Project


# bring in the inp file and convert it to RMark format 
ms_data <- convert.inp("mark_datafiles//do_mark.inp", group.df=data.frame(sex = c("male","female","unidsex")))  #FIXME
                      covariates = data.frame(mass = "sd_mass", guild = c("hgran", "cgran", "foli")))

#---------------------------------------------------------------------------------
#          process multistrata data, includes capture at home, and dipsersal transitions 
#---------------------------------------------------------------------------------
# Build up the model. Looking at sex effects on dispersal/survival
ms_process <- process.data(ms_data,model = "Multistrata",begin.time = 2000,
                           group=c("sex"), sex.var = 1, sex = c(0,1)) #FIXME

ms_ddl <- make.design.data(ms_process)

#---------------------------------------------------------------------------------
#          make dummy variables and covariates
#---------------------------------------------------------------------------------

# Add dynamic dummy variable age class fields to the design data for Psi (transition prob) and p (capture prob)
ms_ddl$Psi$hy = 0
ms_ddl$Psi$hy[ms_ddl$Psi$sex == 0 & ms_ddl$Psi$stratum == "1"] = 1
ms_ddl$Psi$ahy = 0
ms_ddl$Psi$ahy[ms_ddl$Psi$sex >= 1 & ms_ddl$Psi$stratum == "1"] = 1

ms_ddl$p$hy = 0
ms_ddl$p$hy[ms_ddl$p$sex == 1 & ms_ddl$p$stratum == "1"] = 1
ms_ddl$p$hy[ms_ddl$p$sex == 1 & ms_ddl$p$stratum == "1"] = 1
ms_ddl$p$sy = 0
ms_ddl$p$sy[ms_ddl$p$sex == 2 & ms_ddl$p$stratum == "1"] = 1
ms_ddl$p$asy = 0
ms_ddl$p$asy[ms_ddl$p$sex >= 3 & ms_ddl$p$stratum == "1"] = 1
ms_ddl$p$ahy = 0
ms_ddl$p$ahy[ms_ddl$p$sex >= 2 & ms_ddl$p$stratum == "2"] = 1

##### Add dummy variables for operating on specific states or transitions
  # A = 1 (home), B = 2 (away)
  # A to B and B to B, is risky
  # A to A and B to A, is less risky (within home, "normal" movements)
ms_ddl$Psi$toB = 0
ms_ddl$Psi$toB[ms_ddl$Psi$stratum == "1" & ms_ddl$Psi$tostratum == "2"] = 1 

ms_ddl$Psi$toA = 0
ms_ddl$Psi$toA[ms_ddl$Psi$stratum == "2"] = 1

ms_ddl$p$strA = 0
ms_ddl$p$strA[ms_ddl$p$stratum == "1"] = 1
ms_ddl$p$strB = 0
ms_ddl$p$strB[ms_ddl$p$stratum == "2"] = 1

## TODO: fix recapture probabilities for unsampled or omitted months!
    # skipped_periods = c(267, 277, 278, 283, 284, 300, 311, 313, 314, 318, 321, 323, 
              #   337, 339, 344, 351)

# Add a field for monthly reporting probabilities 
ms_ddl$p$rpt=0
ms_ddl$p$rpt[ms_ddl$p$prd] == 267 = 0
ms_ddl$p$rpt[ms_ddl$p$prd] == 277 = 0

# TODO: Lots of work on building up the models!
                          