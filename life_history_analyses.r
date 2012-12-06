# Get life history data on rodent species at Portal

##----- import packages
library(plyr)

##----- functions

##----- set working directory
wd = "C://Users//sarah//Desktop//Dropbox//Active Research Projects//Rodent Movement"
setwd(wd)

# import data
rats = read.csv("PortalRodents_00-09.csv")\

sppList = unique(rats$species)

##----- temporal dynamics

#plot occupancy through time (proportion of years by proportion of months)

for (s in 1:length(sppList)){
  
}


#plot the number of individuals captured in each period across the study

for (s in 1:length(sppList)){
  
}


##----- get fecundity information by species

for (s in 1:length(sppList)){
  sppData = subset(rats, species == sppList[s]) # subset the data by species
  
  females = subset(sppData, sex == "F")
  all_indivs = count(females, vars = "period") #count the number of female individuals in each period
  
  preggers = subset(females, pregnant == "P") #subset data from reproductive females
  repro_indivs = count(preggers, vars = "period")
  
  #join the data together
  d = join(all_indivs, repro_indivs, by = "period", type = "left", match = "all")
  d[is.na(d)] <- 0 #coerce NA to 0
  
  #add column for the proportion of species in each month to be reproductive
  proportion = as.numeric()
  for (row in 1:nrow(d)){
    if (d[row,3] == 0){
      proportion = append(proportion, 0)
    }
    else if (d[row,3] > 0){
      proportion = append(proportion, round(d[row,3]/d[row,2], 3))
    }
  }
  d = cbind(d, proportion) #merge the data
  
  #plot the data
  plot(d$period, d$proportion, xlab = "period", ylab = paste("proportion reproductive", sppList[s], sep = " "), 
       type = "b", pch = 19)
}

