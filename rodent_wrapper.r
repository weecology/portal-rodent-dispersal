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

# give untagged individuals a unique 7-number code
all2 = id_unknowns(all2, 16)

# get rid of 'bad data'; deletes data that is not a pit tag, where sex is inconsistent or where species is inconsistent. 
all2 = subsetDat(all2)

#---------------------------------------------------------------------------------
#          calculate life-history details - reproduction, temporal persistence
#---------------------------------------------------------------------------------

spp_list = c("DO", "DM", "DS", "PI", "PP", "PB", "PF",
             "RM", "RF","RO", "PE", "PP", "PM", "PL", BA",
             "NAO", "SH", "SF", "SO",
             "OT", "OL")

for (spp in 1:length(spp_list)){

}

# average proportion of reproductive females by month across all years

# proportion of reproductive females by month and year

# track the number of times females uniquely reproduce within years

# proportion of years they were seen in

# average number of months they were seen in during years in which they were present

# mean abundance within all years 

# mean abundances within all years on control plots only


#---------------------------------------------------------------------------------
#          calculate movement distances, multi-state capture histories
#---------------------------------------------------------------------------------

# get a vector unique tags, then get a vector of distances moved for all recaptured individuals, by SPECIES
