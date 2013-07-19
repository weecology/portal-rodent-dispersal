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

spplist = c("DO", "DM", "DS", "PI", "PP", "PB", "PF",
             "RM", "RF","RO", "PE", "PM", "PL", "BA",
             "NAO", "SH", "SF",
             "OT", "OL")

for (i in 1:length(spplist)){
  #subset species data
  spdata = subset(all2, species == spplist[i])
  
    # proportion of years they were seen in
    propyrs = length(unique(spdata$yr))/10
      assign(paste(spplist[i], 'propyrs', sep = ""), propyrs)
    # average number of months they were seen in during years in which they were present
    avgmos = mean_win_yr_occ(spdata)
      assign(paste(spplist[i], 'avgmos', sep = ""), avgmos)
    # mean abundance within all years 
    avgabun = allyrs_abun(spdata)
      assign(paste(spplist[i], 'avgabun', sep = ""), avgabun)
    #mean abundances within all years on control plots only
    conabun = allyrs_abun(subset(spdata, plot %in% c(1,2,4,8,9,11,12,14,17,22)))
      assign(paste(spplist[i], 'conabun', sep = ""), conabun)
  
  #subset females for each species
  spdataF = subset(spdata, sex == "F")
  
    #average proportion of reproductive females by month across all years
    reprd = mean_mo_repro(spdataF) #vector with 12 items
      assign(paste(spplist[i], 'propyrs', sep = ""), reprd)
    # proportion of reproductive females by month and year
    reprdyr = mo_repro(spdataF) #matrix with 120 x 5
      assign(paste(spplist[i], 'propyrs', sep = ""), reprodyr)
    # track the number of times individual females uniquely reproduce within years
    indreprd = indiv_repro(spdataF) #matix with cols "tag", "year" and "num_reprod"
      assign(paste(spplist[i], 'propyrs', sep = ""), indreprd)
}


#---------------------------------------------------------------------------------
#          calculate movement distances, multi-state capture histories
#---------------------------------------------------------------------------------

for (i in 1:length(spplist)){
  #subset species data
  spdata = subset(all2, species == spplist[i])
  
    # get a vector unique tags, then get a vector of distances moved for all recaptured individuals, by SPECIES
    tags = unique(spdata$tag)
      assign(paste(spplist[i], 'propyrs', sep = ""), tags)
    meters = distance_moved(spdata, tags)
     assign(paste(spplist[i], 'propyrs', sep = ""), meters)
}

# concatenate core heteromyid data - used to ask if these species behave differently from others
corehet = c(dmmeters, dometers, pbmeters, ppmeters)

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


  