# Code for working with individual-level rodent data
# movement with stakes

library(calibrate)
library(fields)

#---------------------------------------------------------------------------------
#          setup - select wd, import data, source code,  file to collect results
#---------------------------------------------------------------------------------
#set working directory
wd = "/Users/sarah/Documents/GitHub/portal-rodent-dispersal"
setwd(wd)
source("movement_fxns.R")

#import all data
all = read.table('rawdata/all_1980-2009.txt', sep = ',', header = TRUE)

#---------------------------------------------------------------------------------
#          clean up the data
#---------------------------------------------------------------------------------

# change some cols from factor to character class
all$tag = as.character(all$tag)
all$species = as.character(all$species)
all$sex = as.character(all$sex)

#subset data for 10 year of analysis
all2 = subset(all, yr > 1999)
#subset data where species are known (e.g., no "unidentified rodents" or genus-only)
all2 = subset(all2, species!="DX" & species!="UR" & species!="RX" & species!="SX" & species!="PX")

# give untagged individuals a unique 7-number code
all2 = id_unknowns(all2, 16)

# get rid of 'bad data'; deletes data that is not a pit tag, where sex is inconsistent or where species is inconsistent. 
all2 = subsetDat(all2)

#---------------------------------------------------------------------------------
#          calculate life-history details - reproduction, temporal persistence
#---------------------------------------------------------------------------------

spplist = c("DO", "DM", "DS", "PP", "PB", "PF","PI", 
            "PE", "PM", "PL", "RM", "RF","RO", "BA",
            "NAO", "SH", "SF",
            "OT", "OL")

for (i in 1:length(spplist)){
  #subset species data
  spdata = subset(all2, species == spplist[i])
  
    # proportion of years they were seen in
    propyrs = length(unique(spdata$yr))/10
      assign(paste0(spplist[i], 'yr'), propyrs)
    # proportion of years they were seen in CONTROL plots
    conspdata = subset(spdata, plot %in% c(1,2,4,8,9,11,12,14,17,22))
    conpropyrs = length(unique(conspdata$yr))/10
    assign(paste0(spplist[i], 'conyr'), conpropyrs)
    # average number of months they were seen in during years in which they were present
    avgmos = mean_win_yr_occ(spdata)
      assign(paste0(spplist[i], 'avgmos'), avgmos)
    # mean abundance within all years 
    avgabun = allyrs_abun(spdata)
      assign(paste0(spplist[i], 'avgabun'), avgabun)
    # average number of months they were seen in during years in which they were present on CONTROLS
    conavgmos = mean_win_yr_occ(subset(spdata, plot %in% c(1,2,4,8,9,11,12,14,17,22)))
      assign(paste0(spplist[i], 'conavgmos'), conavgmos)
   # number of unique individuals in each year on control plots only
    conabun = allyrs_abun(subset(spdata, plot %in% c(1,2,4,8,9,11,12,14,17,22)))
      assign(paste0(spplist[i], 'conabun'), conabun)
    # mean number of unique individuals on control plots only
    meanconabun = mean(conabun)
      assign(paste0(spplist[i],'meanconabun'), meanconabun)
    # max number of unique individuals on control plots only
    maxconabun = max(conabun)
      assign(paste0(spplist[i], 'maxconabun'), maxconabun)
  
  #subset females for each species
  spdataF = subset(spdata, sex == "F")
  
    #average proportion of reproductive females by month across all years
    reprd = mean_mo_repro(spdataF) #vector with 12 items
      assign(paste0(spplist[i], 'reprd'), reprd)
    # proportion of reproductive females by month and year
    reprdyr = mo_repro(spdataF) #matrix with 120 x 5
      assign(paste0(spplist[i], 'reprdyr'), reprdyr)
    # track the number of times individual females uniquely reproduce within each year
    irep = indiv_repro(spdataF) #matix with cols "tag", "year" and "num_reprod"
      assign(paste0(spplist[i], 'irep'), irep)
}

control_abuns = ls(pattern = "*conabun")
abuns = cbind(BAconabun, DMconabun, DOconabun, DSconabun, NAOconabun, OLconabun, OTconabun, 
              PBconabun, PEconabun, PFconabun, PIconabun, PLconabun, PMconabun, PPconabun,
              RFconabun, RMconabun, ROconabun, SFconabun, SHconabun)
meanabuns = c(BAmeanconabun, DMmeanconabun, DOmeanconabun, DSmeanconabun, NAOmeanconabun, OLmeanconabun, OTmeanconabun,
              PBmeanconabun, PEmeanconabun, PFmeanconabun, PImeanconabun, PLmeanconabun, PMmeanconabun, PPmeanconabun,
              RFmeanconabun, RMmeanconabun, ROmeanconabun, SFmeanconabun, SHmeanconabun)
maxabuns = c(BAmaxconabun, DMmaxconabun, DOmaxconabun, DSmaxconabun, NAOmaxconabun, OLmaxconabun, OTmaxconabun, 
             PBmaxconabun, PEmaxconabun, PFmaxconabun, PImaxconabun, PLmaxconabun, PMmaxconabun, PPmaxconabun, 
             RFmaxconabun, RMmaxconabun, ROmaxconabun, SFmaxconabun, SHmaxconabun)
conyears = cbind(BAconyr, DMconyr, DOconyr, DSconyr, NAOconyr, OLconyr, OTconyr,
                 PBconyr, PEconyr, PFconyr, PIconyr, PLconyr, PMconyr, PPconyr,
                 RFconyr, RMconyr, ROconyr, SFconyr, SHconyr)
conmos = cbind(BAconavgmos, DMconavgmos, DOconavgmos, DSconavgmos, NAOconavgmos, OLconavgmos, OTconavgmos, 
               PBconavgmos, PEconavgmos, PFconavgmos, PIconavgmos, PLconavgmos, PMconavgmos, PPconavgmos, 
               RFconavgmos, RMconavgmos, ROconavgmos, SFconavgmos, SHconavgmos)
#---------------------------------------------------------------------------------
#          calculate movement distances, multi-state capture histories
#---------------------------------------------------------------------------------

for (i in 1:length(spplist)){
  #subset species data
  spdata = subset(all2, species == spplist[i])
  
    # get a vector unique tags, then get a vector of distances moved for all recaptured individuals, by SPECIES
    tags = unique(spdata$tag)
      assign(paste0(spplist[i], 'tags'), tags)
    meters = distance_moved(spdata, tags)
     assign(paste0(spplist[i], 'meters'), meters)
}

# concatenate core heteromyid data - used to ask if these species behave differently from others
corehet = c(DMmeters, DOmeters, PBmeters, PPmeters)

  # find breakpoints to use in MARK data structure for future analyses
  # data reasonably well fits a lognormal distribution (eyeball and J. Powell)
  # breakpoint = mean(logdata) + sd(logdata) of all the distances traveled by recaptured individuals    
  # using log1p, and back transforming using expm1 should solve the problem of lots of zeros 
  corehet_brkpt = expm1(mean(log1p(corehet)) + sd(log1p(corehet)))

# Get MARK capture histories
#------------------------------
periods = c(261:380)
exclosures = c(5, 7, 10, 16, 23, 24)
krat_excl = c(5, 7, 10, 16, 23, 24, 3, 6, 13, 15, 18, 19, 20, 21)

for (i in 1:length(spplist)){
  #subset species data
  spdata = subset(all2, species == spplist[i])
  
  if (i == 1) {
  MARK = noplacelikehome(spdata, periods, krat_excl, corehet_brkpt)
  }
  
  else {
  nextMARK = noplacelikehome(spdata, periods, krat_excl, corehet_brkpt)
  MARK = rbind(MARK, nextMARK)
  }
}

write.table(MARK, file = "mark_datafiles//all_mark.inp", row.names = F, col.names = F, quote = F)


#---------------------------------------------------------------------------------
#          plot results
#---------------------------------------------------------------------------------

#----------------------------------------- plot abundance vs. years, ala core v. transient literature
plot(conyears, maxabuns, pch = 19, xlab = "Persistence (proportion of years present)",
     ylab = "Maximum abundance in any year", bty = "n")
  textxy(conyears, maxabuns, c("BA", "DM", "DO", "DS", "NA", "OL", "OT", "PB", "PE", "PF", "PI",
                              "PL", "PM", "PP", "RF", "RM", "RO", "SF", "SH"), cx = 0.5)

plot(conyears, meanabuns, pch = 19, xlab = "Persistence (proportion of years present)",
     ylab = "Mean yearly abundance across time-series", bty = "n")
textxy(conyears, meanabuns, c("BA", "DM", "DO", "DS", "NA", "OL", "OT", "PB", "PE", "PF", "PI",
                             "PL", "PM", "PP", "RF", "RM", "RO", "SF", "SH"), cx = 0.5)

plot(conyears, conmos, pch = 19, xlab = "Persistence (proportion of years present)", 
     ylab = "Average proportion of months present", xlim = c(0,1), ylim = c(0,1), bty = "n",
     cex.axis = 1.5, cex.lab = 1.5, cex = 1.5)
#      textxy(conyears, conmos, c("BA", "DM", "DO", "DS", "NA", "OL", "OT", "PB", "PE", "PF", "PI",
#                                    "PL", "PM", "PP", "RF", "RM", "RO", "SF", "SH"), cx = 0.5)
    abline(h = 0.66, lty = 2, col = 'gray40', lwd = 1, cex = 1.5)
    abline(v = 0.66, lty = 2, col = 'gray40', lwd = 1, cex = 1.5)


#------------------------------------------ FIGURE - for ESA talk

#plot the relative abundance of females who represent each number of reproductive events per year
plot(NA, NA, xlim = c(0, 4), xaxp = c(0, 4, 4), xlab = "number reproductive events per year", 
     ylab = "proportion females", type = "b", pch = 19, ylim = c(0, 1), yaxp = c(0, 1, 4), bty = "n", 
     main = "core granivores", cex.axis = 1.5, cex.lab = 1.5)
points(c(0:3), table(DOirep$num_reprod)/sum(table(DOirep$num_reprod)), type = "b", pch = 19, col = "black")
points(c(0:4), table(DMirep$num_reprod)/sum(table(DMirep$num_reprod)), type = "b", pch = 19, col = "black")
points(c(0:4), table(PBirep$num_reprod)/sum(table(PBirep$num_reprod)), type = "b", pch = 19, col = "black")
points(c(0:4), table(PPirep$num_reprod)/sum(table(PPirep$num_reprod)), type = "b", pch = 29, col = "black")

plot(NA, NA, xlim = c(0, 4), xaxp = c(0, 4, 4), xlab = "number reproductive events per year", 
     ylab = "proportion females", type = "b", pch = 19, ylim = c(0, 1), yaxp = c(0, 1, 4), bty = "n", 
     main = "transient granivores", cex.axis = 1.5, cex.lab = 1.5)
points(c(0:2), table(PMirep$num_reprod)/sum(table(PMirep$num_reprod)), type = "b", pch = 15, col = "black")
points(c(0:2), table(PEirep$num_reprod)/sum(table(PEirep$num_reprod)), type = "b", pch = 15, col = "black")
points(c(0:1), table(PLirep$num_reprod)/sum(table(PLirep$num_reprod)), type = "b", pch = 15, col = "black")
points(c(0:2), table(PFirep$num_reprod)/sum(table(PFirep$num_reprod)), type = "b", pch = 15, col = "black")
points(c(0:1), table(RMirep$num_reprod)/sum(table(RMirep$num_reprod)), type = "b", pch = 15, col = "black")
points(c(0:1), table(ROirep$num_reprod)/sum(table(ROirep$num_reprod)), type = "b", pch = 15, col = "black")
points(0, table(RFirep$num_reprod)/sum(table(RFirep$num_reprod)), type = "b", pch = 15, col = "black")
points(0, table(DSirep$num_reprod)/sum(table(DSirep$num_reprod)), type = "b", pch = 15, col = "black")
points(c(0:1), table(PIirep$num_reprod)/sum(table(PIirep$num_reprod)), type = "b", pch = 15, col = "black")


#------------------------------------------ FIGURE - for ESA talk

#granivore data from dissertation chapter table 2-3
Psi = c(0.09, 0.13, 0.08, 0.14, 0.48, 0.23, 0.49, 0.41)
S = c(0.76, 0.78, 0.80, 0.81, 0.67, 0.81, 0.76, 0.62)
distbench = c(28.36, 32.64, 25.57, 36.22, 93.05, 74.29, 29.12, 93.30)

# S vs. Psi
plot(Psi, S, pch = 19, xlim = c(0:1), ylim = c(0:1), bty = "n", cex.axis = 1.5, cex.lab = 1.5, cex = 1.5,
     ylab = "survival probability", xlab = "long distance movement probability")
plot(distbench, S, pch = 19, bty = "n", xlab = "transition distance (meters)", ylab = "survival", 
     xlim = c(0,100), ylim = c(0,1), cex.axis = 1.5, cex.lab = 1.5, cex = 1.5)

lm1 = lm(S~Psi)
lm2 = lm(S~distbench)


  