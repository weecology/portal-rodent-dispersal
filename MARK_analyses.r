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
SA = as.numeric(row.names(MS.ddl$S[MS.ddl$S$stratum == "A",]))
SB = as.numeric(row.names(MS.ddl$S[MS.ddl$S$stratum == "B",]))
SC = as.numeric(row.names(MS.ddl$S[MS.ddl$S$stratum == "C",]))  #No stratum C??
SAval = rep(1,length(SA))
SBval = rep(0,length(SB))
SCval = rep(0,length(SC))
# fix S to 1 for A and 0 for dead states B and C
Sfix <- list(formula=~stratum,fixed=list(index=c(SA,SB,SC),value=c(SAval,SBval,
                                                                   SCval)))

#---------------------------------------------------------------------------------
#          Define model structures for p (capture probability)
#---------------------------------------------------------------------------------

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




