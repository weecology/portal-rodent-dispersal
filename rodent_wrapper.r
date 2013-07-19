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
    # average number of months they were seen in during years in which they were present
    avgmos = mean_win_yr_occ(spdata)
    # mean abundance within all years 
    avgabun = allyrs_abun(spdata)
    #mean abundances within all years on control plots only
    conabun = allyrs_abun(subset(spdata, plot %in% c(1,2,4,8,9,11,12,14,17,22)))
  
  #subset females for each species
  spdataF = subset(spdata, sex == "F")
  
    #average proportion of reproductive females by month across all years
    reprd = mean_mo_repro(spdataF) #vector with 12 items
    # proportion of reproductive females by month and year
    reprdyr = mo_repro(spdataF) #matrix with 120 x 5
    # track the number of times individual females uniquely reproduce within years
    indreprd = indiv_repro(spdataF) #matix with cols "tag", "year" and "num_reprod"

}


# average number of months they were seen in during years in which they were present

# mean abundance within all years 

# mean abundances within all years on control plots only


#---------------------------------------------------------------------------------
#          calculate movement distances, multi-state capture histories
#---------------------------------------------------------------------------------

# get a vector unique tags, then get a vector of distances moved for all recaptured individuals, by SPECIES
