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


# bring in the inp file and convert it to RMark format - This file includes all the data from all the species 
ms_data <- convert.inp("mark_datafiles//all_mark.inp", group.df=data.frame(sex = c("male","female","unidsex")),  
                      covariates = c("sd_mass", "guild", "species"))

#---------------------------------------------------------------------------------
#          process multistrata data, includes capture at home, and dipsersal transitions 
#---------------------------------------------------------------------------------
# Build up the model. Looking at sex effects on dispersal/survival
ms_process <- process.data(ms_data,model = "Multistrata", begin.time = 261, groups = "sex")

ms_ddl <- make.design.data(ms_process)

#---------------------------------------------------------------------------------
#          make dummy variables and covariates
#---------------------------------------------------------------------------------
# Add dummy variables for operating on specific states or transitions
# A = 1 (home), B = 2 (away)
# A to B and B to B, is risky
# A to A and B to A, is less risky (within home, "normal" movements)

# Surival probability given that the individual is in A
ms_ddl$S$inA = 0
ms_ddl$S$inA[ms_ddl$S$stratum == "1"] = 1

# Surival probability given that the individual is in B
ms_ddl$S$inB = 0
ms_ddl$S$inB[ms_ddl$S$stratum == "2"] = 1

# Transition probability given that the individual  A ---> B
ms_ddl$Psi$toA = 0
ms_ddl$Psi$toA[ms_ddl$Psi$stratum == "2" & ms_ddl$Psi$tostratum == "1"] = 1

# Transition probability given that the individual  B ---> A
ms_ddl$Psi$toB = 0
ms_ddl$Psi$toB[ms_ddl$Psi$stratum == "1" & ms_ddl$Psi$tostratum == "2"] = 1 

# recapture probability given that the individual is in A
ms_ddl$p$strA = 0
ms_ddl$p$strA[ms_ddl$p$stratum == "1"] = 1

# recapture probability given that the individual is in B
ms_ddl$p$strB = 0
ms_ddl$p$strB[ms_ddl$p$stratum == "2"] = 1


#--------------------------------------------------------------------------------
    # TODO: Build up the models

    # Do psi and S significantly differ among species/guilds?

    # Within a species/guild, do covariates (sex, body mass) influence psi and S?

#---------------------------------------------------------------------------------
#          Define model structures for S (survival probability)
#---------------------------------------------------------------------------------
Snull <- list(formula = ~1)           # null model, S is not dependent on strata
Sstrata <- list(formula = ~stratum)   # S is dependent on strata (in A or in B)

Sguild <- list(formula = ~guild)  # Is this how to set up a model testing if S differs among guilds?
Sspecies <- list(formula = ~species) # test for differences among species

#---------------------------------------------------------------------------------
#          Define model structures for p (capture probability)
#---------------------------------------------------------------------------------
# fix recapture probabilities for unsampled or omitted months
#    skipped_periods = c(267, 277, 278, 283, 284, 300, 311, 313, 314, 318, 321, 323, 337, 339, 344, 351): p = 0

# select periods that were omitted from the study - untrapped  #FIXME, error in undefined columns selected?
p267 = as.numeric(row.names(ms_ddl$p[ms_ddl$p$time == 267]))
p277 = as.numeric(row.names(ms_ddl$p[ms_ddl$p$time == 277]))
p278 = as.numeric(row.names(ms_ddl$p[ms_ddl$p$time == 278]))
p283 = as.numeric(row.names(ms_ddl$p[ms_ddl$p$time == 283]))
p284 = as.numeric(row.names(ms_ddl$p[ms_ddl$p$time == 284]))
p300 = as.numeric(row.names(ms_ddl$p[ms_ddl$p$time == 300]))
p311 = as.numeric(row.names(ms_ddl$p[ms_ddl$p$time == 311]))
p313 = as.numeric(row.names(ms_ddl$p[ms_ddl$p$time == 313]))
p314 = as.numeric(row.names(ms_ddl$p[ms_ddl$p$time == 314]))
p318 = as.numeric(row.names(ms_ddl$p[ms_ddl$p$time == 318]))
p321 = as.numeric(row.names(ms_ddl$p[ms_ddl$p$time == 321]))
p323 = as.numeric(row.names(ms_ddl$p[ms_ddl$p$time == 323]))
p337 = as.numeric(row.names(ms_ddl$p[ms_ddl$p$time == 337]))
p339 = as.numeric(row.names(ms_ddl$p[ms_ddl$p$time == 339]))
p344 = as.numeric(row.names(ms_ddl$p[ms_ddl$p$time == 344]))
p351 = as.numeric(row.names(ms_ddl$p[ms_ddl$p$time == 351]))

# set those periods to p = 0, because they *can't* be anything else
p267val = rep(0, length(p267))
p277val = rep(0, length(p277))
p278val = rep(0, length(p278))
p283val = rep(0, length(p283))
p284val = rep(0, length(p284))
p300val = rep(0, length(p300))
p311val = rep(0, length(p311))
p313val = rep(0, length(p313))
p314val = rep(0, length(p314))
p318val = rep(0, length(p318))
p321val = rep(0, length(p321))
p323val = rep(0, length(p323))
p337val = rep(0, length(p337))
p339val = rep(0, length(p339))
p344val = rep(0, length(p344))
p351val = rep(0, length(p351))


# look for effects of strata on recapture probability, given that some p are fixed to 0 (listed below)
pstrata <- list(formula = ~stratum, fixed = list(index = c(p267, p277, p278, p283, p284, p300, p311, p313, p314,
                                                     p318, p321, p323, p337, p339, p344, p351), 
                                             value = c(p267val, p277val, p278val, p283val, p284val, p300val, p311val,
                                                       p313val, p314val, p318val, p321val, p323val, p337val, p339val,
                                                       p344val, p351val), link = "cloglog"))
        # link = "logit" is the default. "cloglog" may be esp. useful when there are fewer recaptures


pguild <- list(formula = ~guild, fixed = list(index = c(p267, p277, p278, p283, p284, p300, p311, p313, p314,
                                                         p318, p321, p323, p337, p339, p344, p351), 
                                               value = c(p267val, p277val, p278val, p283val, p284val, p300val, p311val,
                                                         p313val, p314val, p318val, p321val, p323val, p337val, p339val,
                                                         p344val, p351val), link = "cloglog"))
              # Is this how to set up a model testing if S differs among guilds?


#---------------------------------------------------------------------------------
#          Define model structures for Psi (transition probability)
#---------------------------------------------------------------------------------
Psistrata <- list(formula = ~stratum)

Psiguild <- list(formula = ~guild)


#---------------------------------------------------------------------------------
#          Run Models and collect results
#---------------------------------------------------------------------------------
#SIMANNEAL should be best for multistrata models, but may take longer to run
Sstrata_pstrata_Psistrata <- mark(ms_process,ms_ddl, model.parameters = list(S = Sstrata,  p = pstrata, Psi = Psistrata),
                                  options = "SIMANNEAL")

Sguild_pguild_Psiguild <- mark(ms_process,ms_ddl, model.parameters = list(S = Sguild,  p = pguild, Psi = Psiguild),
                               options = "SIMANNEAL")

Sguild_pstrata_Psistrata <- mark(ms_process,ms_ddl, model.parameters = list(S = Sguild,  p = pstrata, Psi = Psiguild),
                                 options = "SIMANNEAL")

Sstrata_pstrata_Psiguild <- mark(ms_process,ms_ddl, model.parameters = list(S = Sstrata,  p = pstrata, Psi = Psiguild),
                                 options = "SIMANNEAL")

Sguild_pstrata_Psistrata <- mark(ms_process,ms_ddl, model.parameters = list(S = Sguild,  p = pstrata, Psi = Psistrata),
                                 options = "SIMANNEAL")

Sspecies_pstrata_Psistrata <- mark(ms_process,ms_ddl, model.parameters = list(S = Sspecies,  p = pstrata, Psi = Psistrata),
                                   options = "SIMANNEAL")


#summarize results
ms_results <- collect.models(type = "Multistrata")


#---------------------------------------------------------------------------------
#          Write result data to csv files
#---------------------------------------------------------------------------------
write.csv(Sstrata_pstrata_Psistrata$results$beta, "ms_modelstuff_beta")
write.csv(Sstrata_pstrata_Psistrata$results$real, "ms_modelstuff_real")
#TODO: Add the others to be printed when finalize models

