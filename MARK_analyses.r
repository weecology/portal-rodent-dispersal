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
ms_process <- process.data(ms_data,model = "Multistrata",begin.time = 261, groups = "sex")

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
    # TODO: Lots of work on building up the models!
        
    # What is psi (transition probability) for each species/guild?

    # What is  S (survival probability) for each species/guild?

    # What is p (capture probability for each species/guild?

    # Do psi and S significantly differ among species/guilds?

    # Within a species/guild, do covariates (sex, body mass) influence psi and S?

#---------------------------------------------------------------------------------
#          Define model structures for S (survival probability)
#---------------------------------------------------------------------------------
Snull <- list(formula=~1)           # null model, S is not dependent on strata
Sstrata <- list(formula=~stratum)   # S is dependent on strata (in A or in B)

#---------------------------------------------------------------------------------
#          Define model structures for p (capture probability)
#---------------------------------------------------------------------------------

#---------------------------------------------------------------------------------
#          fix p for omitted  periods - time between  trapping events ~ 1 month
#---------------------------------------------------------------------------------
## TODO: fix recapture probabilities for unsampled or omitted months!
# skipped_periods = c(267, 277, 278, 283, 284, 300, 311, 313, 314, 318, 321, 323, 
#   337, 339, 344, 351): p = 0

# Add a field for monthly reporting probabilities 
ms_ddl$p$rpt = 0
ms_ddl$p$rpt[ms_ddl$p$prd] == 267 = 0  #FIX ME: is this on the right track?
ms_ddl$p$rpt[ms_ddl$p$prd] == 277 = 0

  ## FIXME: What does this all mean??
p1996=as.numeric(row.names(MS.ddl$p[MS.ddl$p$time==1996&MS.ddl$p$stratum=="A",]))
p1997=as.numeric(row.names(MS.ddl$p[MS.ddl$p$time==1997&MS.ddl$p$stratum=="A",]))
p1998.g=as.numeric(row.names(MS.ddl$p[MS.ddl$p$time==1998&MS.ddl$p$stratum=="A"&MS.ddl$p$rage=="gosling",]))
p1998.y=as.numeric(row.names(MS.ddl$p[MS.ddl$p$time==1998&MS.ddl$p$stratum=="A"&MS.ddl$p$rage=="yearling",]))
p2003.g=as.numeric(row.names(MS.ddl$p[MS.ddl$p$time==2003&MS.ddl$p$stratum=="A"&MS.ddl$p$rage=="gosling",]))
p2005.g=as.numeric(row.names(MS.ddl$p[MS.ddl$p$time==2005&MS.ddl$p$stratum=="A"&MS.ddl$p$rage=="gosling",]))
p2009=as.numeric(row.names(MS.ddl$p[MS.ddl$p$time==2009&MS.ddl$p$stratum=="A",]))
p2010.g=as.numeric(row.names(MS.ddl$p[MS.ddl$p$time==2010&MS.ddl$p$stratum=="A"&MS.ddl$p$rage=="gosling",]))
rC=as.numeric(row.names(MS.ddl$p[MS.ddl$p$stratum=="C",]))
r1997.g=as.numeric(row.names(MS.ddl$p[MS.ddl$p$time==1997&MS.ddl$p$stratum=="B"&MS.ddl$p$rage=="gosling",]))
r1998.g=as.numeric(row.names(MS.ddl$p[MS.ddl$p$time==1998&MS.ddl$p$stratum=="B"&MS.ddl$p$rage=="gosling",]))
r2003.g=as.numeric(row.names(MS.ddl$p[MS.ddl$p$time==2003&MS.ddl$p$stratum=="B"&MS.ddl$p$rage=="gosling",]))
r2005.g=as.numeric(row.names(MS.ddl$p[MS.ddl$p$time==2005&MS.ddl$p$stratum=="B"&MS.ddl$p$rage=="gosling",]))
r2010.g=as.numeric(row.names(MS.ddl$p[MS.ddl$p$time==2010&MS.ddl$p$stratum=="B"&MS.ddl$p$rage=="gosling",]))
p1996val=rep(0,length(p1996))
p1997val=rep(0,length(p1997))
p1998.gval=rep(0,length(p1998.g))
p1998.yval=rep(0,length(p1998.y))
p2003.gval=rep(0,length(p2003.g))
p2005.gval=rep(0,length(p2005.g))
p2009val=rep(0,length(p2009))
p2010.gval=rep(0,length(p2010.g))
rCval=rep(0,length(rC))
r1997.gval=rep(0,length(r1997.g))
r1998.gval=rep(0,length(r1998.g))
r2003.gval=rep(0,length(r2003.g))
r2005.gval=rep(0,length(r2005.g))
r2010.gval=rep(0,length(r2010.g))

#---------------------------------------------------------------------------------
#          Define model structures for psi (transition probability)
#---------------------------------------------------------------------------------
## FIXME: What does this all mean??
PsiB = as.numeric(row.names(MS.ddl$Psi[MS.ddl$Psi$stratum == "B",]))
PsiC = as.numeric(row.names(MS.ddl$Psi[MS.ddl$Psi$stratum == "C",]))
Psi1996.g = as.numeric(row.names(MS.ddl$Psi[MS.ddl$Psi$time == 1996&MS.ddl$Psi$stratum == "A"&MS.ddl$Psi$rage == "gosling",]))
Psi1997.g = as.numeric(row.names(MS.ddl$Psi[MS.ddl$Psi$time == 1997&MS.ddl$Psi$stratum == "A"&MS.ddl$Psi$rage == "gosling",]))
PsiBval = rep(0,length(PsiB))
PsiCval = rep(0,length(PsiC))
Psi1996.gval = rep(0,length(Psi1996.g))
Psi1997.gval = rep(0,length(Psi1997.g))
Psi2002.gval = rep(0,length(Psi2002.g))
Psi2004.gval = rep(0,length(Psi2004.g))
Psi2009.gval = rep(0,length(Psi2009.g))

# tostratum and age effects (unique age effects for harvest and non-harvest mortalities)
Psiage <- list(formula=~-1+toB:hy+toB:ahy+toC:hy+toC:ahy+fromB+fromC,
               fixed=list(index=c(PsiB,PsiC,Psi1996.g,Psi1997.g,Psi2002.g,Psi2004.g,Psi2009.g),
                          value=c(PsiBval,PsiCval,Psi1996.gval,Psi1997.gval,Psi2002.gval,Psi2004.gval,
                                  Psi2009.gval)),link="logit")

#---------------------------------------------------------------------------------
#          Run Models
#---------------------------------------------------------------------------------
Sfix.prage.Psiage <- mark(MS.process,MS.ddl,model.parameters=list(S=Sfix,
                                                                  p=prage,Psi=Psiage),options="SIMANNEAL")
Sfix.ptbin.Psiage <- mark(MS.process,MS.ddl,model.parameters=list(S=Sfix,
                                                                  p=ptbin,Psi=Psiage),options="SIMANNEAL")
Sfix.prtbin.Psiage <- mark(MS.process,MS.ddl,model.parameters=list(S=Sfix,
                                                                   p=prtbin,Psi=Psiage),options="SIMANNEAL")
Sfix.ptbinrT3.Psiage <- mark(MS.process,MS.ddl,model.parameters=list(S=Sfix,
                                                                     p=ptbinrT3,Psi=Psiage))
Sfix.ptbinrrpt.Psiage <- mark(MS.process,MS.ddl,model.parameters=list(S=Sfix,
                                                                      p=ptbinrrpt,Psi=Psiage),options="SIMANNEAL")
Sfix.pT3.Psiage <- mark(MS.process,MS.ddl,model.parameters=list(S=Sfix,
                                                                p=pT3,Psi=Psiage),options="SIMANNEAL")
Sfix.pT3rtbin.Psiage <- mark(MS.process,MS.ddl,model.parameters=list(S=Sfix,
                                                                     p=pT3rtbin,Psi=Psiage))
Sfix.prT3.Psiage <- mark(MS.process,MS.ddl,model.parameters=list(S=Sfix,
                                                                 p=prT3,Psi=Psiage),options="SIMANNEAL")
Sfix.pT3rrpt.Psiage <- mark(MS.process,MS.ddl,model.parameters=list(S=Sfix,
                                                                    p=pT3rrpt,Psi=Psiage),options="SIMANNEAL")

#summarize results
MS.results <- collect.models(type="Multistrata")

#---------------------------------------------------------------------------------
#          Write result data to csv files
#---------------------------------------------------------------------------------
write.csv(Sfix.pT3rrpt.Psiagebytbin$results$beta, "MS_agebytbin_beta")
write.csv(Sfix.pT3rrpt.Psiagebytbin$results$real, "MS_agebytbin_real")
write.csv(Sfix.pT3rrpt.Psiagebyt$results$beta, "MS_agebyt_beta")
write.csv(Sfix.pT3rrpt.Psiagebyt$results$real, "MS_agebyt_real")