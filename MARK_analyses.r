# for manipulating  rodent movement data in MARK
## analyze survival and dispseral probabilities for species and guilds
# Use multistrata models that give S(survival), p(capture probability) and Psi(transition probability)
# Will compare data among species and guilds. 
# Major covariates include sex, guild, and average body mass

# Do psi and S significantly differ among species/guilds?

# Within a species/guild, do covariates (sex, body mass) influence psi and S?

rm(list=ls(all=TRUE))   # clears the memory

#---------------------------------------------------------------------------------
#          bring in the data and source files
#---------------------------------------------------------------------------------
#set working directory and import source code
#setwd("~/")
setwd("~/portal-rodent-dispersal/")


# Run the line below to generate new .inp files 
#source("stake_movement.r") #makes a mark data structure using species-level data from Portal Project

#grab all the .inp files to loop over for analysis
files = list.files(getwd(), pattern = "mark.inp", full.name=T, recursive=T)
files = files[c(5,6,10)]


for (f in 1:length(files)){
  
  require(RMark)
  
# bring in the inp file and convert it to RMark format - This file includes all the data from all the species 
# files are carn_mark.inp, foli_mark.inp and gran_mark.inp
ms_data <- convert.inp(files[f], group.df = data.frame(sex=c("male", "female", "unidentified")),
                      covariates = c("species"))

#convert to factor
#ms_data$guild = as.factor(ms_data$guild)
ms_data$species = as.factor(ms_data$species)
#ms_data$status = as.factor(ms_data$status)

cat("Imported data.", file="outfile.txt", sep="\n")

#---------------------------------------------------------------------------------
#          process multistrata data, includes capture at home, and dipsersal transitions 
#---------------------------------------------------------------------------------
# Build up the model. Looking at sex effects on dispersal/survival
# begin.time == first period number
ms_process <- process.data(ms_data, model = "Multistrata", begin.time = 130, 
                           groups = c("sex", "species"))

ms_ddl <- make.design.data(ms_process) #ddl = design data list

#---------------------------------------------------------------------------------
#          make dummy variables and covariates
#---------------------------------------------------------------------------------
# Add dummy variables for operating on specific states(strata) or transitions
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
ms_ddl$p$strataA = 0
ms_ddl$p$strataA[ms_ddl$p$stratum == "1"] = 1

# recapture probability given that the individual is in B
ms_ddl$p$strataB = 0
ms_ddl$p$strataB[ms_ddl$p$stratum == "2"] = 1


#--------------------------------------------------------------------------------
#           Build up the models
#---------------------------------------------------------------------------------
#          Define model structures for S (survival probability)
#---------------------------------------------------------------------------------
Snull <- list(formula = ~1)          


#---------------------------------------------------------------------------------
#          Define model structures for p (capture probability)
#---------------------------------------------------------------------------------
# fix recapture probabilities for unsampled or omitted months
#    skipped_periods = c(237, 241, 267, 277, 278, 283, 284, 300, 311, 313, 314, 318, 321, 323, 337, 339, 344, 351): p = 0

# select periods that were omitted from the study - untrapped
p237 = as.numeric(row.names(ms_ddl$p[ms_ddl$p$time == 237,]))
p241 = as.numeric(row.names(ms_ddl$p[ms_ddl$p$time == 241,]))
p267 = as.numeric(row.names(ms_ddl$p[ms_ddl$p$time == 267,])) 
p277 = as.numeric(row.names(ms_ddl$p[ms_ddl$p$time == 277,]))
p278 = as.numeric(row.names(ms_ddl$p[ms_ddl$p$time == 278,]))
p283 = as.numeric(row.names(ms_ddl$p[ms_ddl$p$time == 283,]))
p284 = as.numeric(row.names(ms_ddl$p[ms_ddl$p$time == 284,]))
p300 = as.numeric(row.names(ms_ddl$p[ms_ddl$p$time == 300,]))
p311 = as.numeric(row.names(ms_ddl$p[ms_ddl$p$time == 311,]))
p313 = as.numeric(row.names(ms_ddl$p[ms_ddl$p$time == 313,]))
p314 = as.numeric(row.names(ms_ddl$p[ms_ddl$p$time == 314,]))
p318 = as.numeric(row.names(ms_ddl$p[ms_ddl$p$time == 318,]))
p321 = as.numeric(row.names(ms_ddl$p[ms_ddl$p$time == 321,]))
p323 = as.numeric(row.names(ms_ddl$p[ms_ddl$p$time == 323,]))
p337 = as.numeric(row.names(ms_ddl$p[ms_ddl$p$time == 337,]))
p339 = as.numeric(row.names(ms_ddl$p[ms_ddl$p$time == 339,]))
p344 = as.numeric(row.names(ms_ddl$p[ms_ddl$p$time == 344,]))
p351 = as.numeric(row.names(ms_ddl$p[ms_ddl$p$time == 351,]))

# set those periods to p = 0, because they *can't* be anything else
p237val = rep(0, length(p237))
p241val = rep(0, length(p241))
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

cat("Fixed omitted periods to zero.", sep="\n", file="outfile.txt", append=TRUE)


# look for effects on recapture probability, given that some p are fixed to 0 (listed below)
# link = "logit" is the default. "cloglog" may be esp. useful when there are fewer recaptures

#Null Model
pnull <- list(formula = ~1, fixed = list(index = c(p237, p241, p267, p277, p278, p283, p284, p300, p311, p313, p314,
                                                   p318, p321, p323, p337, p339, p344, p351), 
                                         value = c(p237, p241, p267val, p277val, p278val, p283val, p284val, p300val, p311val,
                                                   p313val, p314val, p318val, p321val, p323val, p337val, p339val,
                                                   p344val, p351val), link = "cloglog"))

cat("Searched for period effect on recapture probability.", sep="\n", file="outfile.txt", append=TRUE)

#---------------------------------------------------------------------------------
#          Define model structures for Psi (transition probability)
#---------------------------------------------------------------------------------
Psinull <- list(formula = ~1, link = "mlogit")

cat("Defined model structure for Psi", sep="\n", file="outfile.txt", append=TRUE)

  
#---------------------------------------------------------------------------------
#          Run Models and collect results
#---------------------------------------------------------------------------------
#send results to new folder - change working directory
wd = "~/portal-rodent-dispersal/mark_output/"
 # wd = "C:/Users/sarah/Documents/GitHub/portal-rodent-dispersal/mark_output"
  setwd(wd)

cat("running the multistrata models", sep="\n", file="outfile.txt", append=TRUE)

# #SIMANNEAL should be best for multistrata models, but may take longer to run
Snull_pnull_Psinull <- mark(ms_process, ms_ddl, model.parameters = list(S=Snull,  p=pnull, Psi=Psinull),
                            options="SIMANNEAL", external=TRUE)

cat("Null model is finished", sep="\n", file="outfile.txt", append=TRUE)


#-----------------------------------------------------
#             summarize results
#-----------------------------------------------------
ms_results <- collect.models(type = "Multistrata")

cat("Summarized results.", sep="\n", file="outfile.txt", append=TRUE)

  print(ms_results)
  print (Snull_pnull_Psinull$results$beta[1:3,])

#---------------------------------------------------------------------------------
#          Write result data to csv files
#---------------------------------------------------------------------------------

write.csv(Snull_pnull_Psinull$results$beta, paste("ms_null_beta_",ms_data[1,6],".csv",sep=""))
write.csv(Snull_pnull_Psinull$results$real, paste("ms_null_real_",ms_data[1,6],".csv",sep=""))

cat("End Code. Look for your csv files.", sep="\n", file="outfile.txt", append=TRUE)
print( paste("file", files[f], " is done.", sep = ""))
  
rm(list=ls()[!ls() %in% c("f", "files")])   # clears the memory of everything except the file list

}
