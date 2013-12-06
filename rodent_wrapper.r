# Code for working with individual-level rodent data
# movement with stakes

library(reshape)
library(ggplot2)
library(calibrate)
library(fields)

#---------------------------------------------------------------------------------
#          setup - select wd, import data, source code,  file to collect results
#---------------------------------------------------------------------------------
#set working directory
wd = "/Users/sarah/Documents/GitHub/portal-rodent-dispersal"
wd = "C:/Users/sarah/Documents/GitHub/portal-rodent-dispersal"
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

allmonths = count_months(all, c(1980:2009))

#subset data where species are known (e.g., no "unidentified rodents" or genus-only)
all2 = subset(all, species!="DX" & species!="UR" & species!="RX" & species!="SX" & species!="PX" & species != "OX")

# give untagged individuals a unique 7-number code
all2 = id_unknowns(all2, 16)

# make sure when note2 == * it is counted as a new tag
tags = unique(all2$tag)
all3 = starred_tags(all2, tags)

# check for individuals with same tag, but captured over long timespan (may be able to separate them) 
# necessary if using ALL data (ear and toe tags)
# returns the dataset with new IDs for those that can easily be sorted out based on sex and year

dups = is_duplicate_tag(all2, tags, 10, 9, 16) #check to see if can separate tags based

#subset data for years of analysis
all2 = subset(all, yr > 1988)

# get rid of 'bad data'; deletes data where species is inconsistent. 
all2 = subsetDat(all2)


#---------------------------------------------------------------------------------
#          calculate life-history details - temporal persistence, reproduction
#---------------------------------------------------------------------------------
spplist = unique(all2$species)
totalyears = c(1989:2009) #will need to be adjusted based on time frame want to use
controlplots = c(1,2,4,8,9,11,12,14,17,22)
months = count_months(all2, totalyears)

persistence = data.frame(species=NA, propyrs=NA, propmos=NA, meanabun=NA, maxabun=NA)
  outcount = 1
yearly_control_abundance = data.frame(year = totalyears)
reprod_mo_yr = data.frame(year=NA, month=NA, proprepro=NA, numfemales=NA, species=NA)
avg_mo_reprod = data.frame(species=NA, month=NA, proprepro=NA)

for (i in 1:length(spplist)){
  #subset species data
  spdata = subset(all2, species == spplist[i])
  
    # proportion of years they were seen in CONTROL plots
    conspdata = subset(spdata, plot %in% controlplots)
    conpropyrs = length(unique(conspdata$yr))/length(totalyears)
    # average proportion of months they were seen in during years in which they were present on CONTROLS
    conavgmos = mean_win_yr_occ(conspdata, totalyears, months)
    # number of unique individuals in each year on control plots only
    conabun = allyrs_abun(conspdata, totalyears)
    # mean number of unique individuals on control plots only, includes ZEROES
    meanconabun = mean(conabun)
    # max number of unique individuals on control plots only
    maxconabun = max(conabun)
  
  #record persistence and abundance data from control plots
  persistence[outcount,] = c(spplist[i], round(conpropyrs,4), round(conavgmos,4), round(meanconabun,4), maxconabun)
  yearly_control_abundance = cbind(yearly_control_abundance, conabun)
  outcount = outcount + 1
  
  #subset females for each species
  spdataF = subset(spdata, sex == "F")
  
    #average proportion of reproductive females by month across all years
    reprd = mean_mo_repro(spdataF, totalyears) #vector with 12 items (one for each month)
      avg_mo_reprod = rbind(avg_mo_reprod, reprd)    
  
  # proportion of reproductive females by month and year
    reprdyr = mo_repro(spdataF) #matrix with 120 x 5
      reprod_mo_yr = rbind(reprod_mo_yr, reprdyr)
  
    # track the number of times individual females uniquely reproduce within each year
    irep = indiv_repro(spdataF) #matix with cols "tag", "year" and "num_reprod"
      assign(paste0(spplist[i], 'irep'), irep)
}

#add species names to the dataframe of abundance vectors
names(yearly_control_abundance) = c("year", spplist)
#make sure cols are numeric
persistence$propyrs = as.numeric(persistence$propyrs)
persistence$propmos = as.numeric(persistence$propmos)
persistence$meanabun = as.numeric(persistence$meanabun)
persistence$maxabun = as.numeric(persistence$maxabun)
avg_mo_reprod$proprepro = as.numeric(avg_mo_reprod$proprepro)
avg_mo_reprod$month = as.numeric(avg_mo_reprod$month)
#delete Null rows
avg_mo_reprod = avg_mo_reprod[-1,]
reprod_mo_yr = reprod_mo_yr[-1,]


#melt control abundance data frame for later plotting
yrcontrolabuns = melt(yearly_control_abundance, id.vars=c("year"))
names(yrcontrolabuns) = c("year", "species", "abun")


#Identify Core species based on results
corespecies = persistence[which(persistence$propyrs >= 0.666 & persistence$propmos >= 0.666),1]
intermediatespecies = persistence[which(persistence$propyrs >= 0.666 & persistence$propmos < 0.666),1]
transientspecies = persistence[which(persistence$propyrs < 0.666),1]


#---------------------------------------------------------------------------------
#          calculate movement distances, multi-state capture histories
#---------------------------------------------------------------------------------

#FIXME: Create an empty list, add each species-level vector as an element in the list
# Then don't have to worry about locating the data later

for (i in 1:length(spplist)){
  #subset species data
  spdata = subset(all2, species == spplist[i])
  
    # get a vector unique tags, then get a vector of distances moved for all recaptured individuals, by SPECIES
    tags = unique(spdata$tag)
      assign(paste0(spplist[i], 'tags'), tags)
    mtrs = distance_moved(spdata, tags)
     assign(paste0(spplist[i], 'meters'), mtrs)
}

# concatenate core granivore data - used to ask if these species behave differently from others
corehet = c(DMmeters, DOmeters, PBmeters, PPmeters)

  # find breakpoints to use in MARK data structure for future analyses
  # data reasonably well fits a lognormal distribution (eyeball and J. Powell)
  # breakpoint = mean(logdata) + sd(logdata) of all the distances traveled by recaptured individuals    
  # using log1p, and back transforming using expm1 should solve the problem of lots of zeros 
  corehet_brkpt = expm1(mean(log1p(corehet)) + sd(log1p(corehet)))

# Get MARK capture histories
#------------------------------
#periods = c(261:380) #2000-2009
periods = c(130:380) #1989-2009
all_excl = c(5, 7, 10, 16, 23, 24) 
krat_excl = c(5, 7, 10, 16, 23, 24, 3, 6, 13, 15, 18, 19, 20, 21)

for (i in 1:length(spplist)){
  #subset species data, for each species in turn
  spdata = subset(all2, species == spplist[i])
  
  if (spplist[i] %in% list("DM", "DS", "DO")) { exclosures = krat_excl} 
  else { exclosures = all_excl} #TOOD: OK unless good reason to only use control plots
  
  #the first species begins the new data matrix for MARK
  if (i == 1) {
  MARK = noplacelikehome(spdata, periods, exclosures, corehet_brkpt)
  }
  
  #all subsequent species are appended onto the end of the existing MARK data matrix
  else {
  nextMARK = noplacelikehome(spdata, periods, exclosures, corehet_brkpt)
  MARK = rbind(MARK, nextMARK)
  }
}

write.table(MARK, file = "mark_datafiles//all_mark.inp", row.names = F, col.names = F, quote = F)


#---------------------------------------------------------------------------------
#          plot results
#---------------------------------------------------------------------------------
#----------------------------- plot abundance vs. years, ala core v. transient literature
ggplot(persistence, aes(propyrs, propmos)) + geom_point(aes(size = maxabun)) + theme_bw() + xlab("proportion years present") +
  ylab("proportion of months present") + xlim(0,1) + ylim(0,1) + 
  geom_vline(xintercept=0.66, linetype="dotted", col = "red") + 
  geom_hline(yintercept=0.66, linetype="dotted", col = "red") +
  ggtitle("Rodents 1989 - 2009")

ggplot(persistence, aes(propyrs, propmos)) + geom_point(aes(size = meanabun)) + theme_bw() + xlab("proportion years present") +
  ylab("proportion of months present") + xlim(0,1) + ylim(0,1) + 
  geom_vline(xintercept=0.66, linetype="dotted", col = "red") + 
  geom_hline(yintercept=0.66, linetype="dotted", col = "red") +
  ggtitle("Rodents 1989 - 2009") + geom_text(aes(label = species), hjust=0, vjust=0)

#------------------------- plot abundance for all species across timeseries
guild = c(3,1,1,1,2,1,2,3,1,1,1,2,1,1,1,1,2,1,1,1,1) #TODO: Color by guild? Add guild to df?
ggplot(yrcontrolabuns, aes(x = year, y = abun, group = species)) + 
  geom_point() + geom_line()

#------------------------- plot monthly reproduction
ggplot(avg_mo_reprod, aes(month, proprepro)) + geom_point() + 
  geom_line() + facet_wrap(~species)

#------------------------- plot meters traveled by all species
distances = ls(pattern = "*meters") #see all the vectors

ggplot(data=data.frame(value=DOmeters)) +
  stat_density(aes(x=value)) +
  geom_segment(aes(x=value,xend=value),y=0,yend=0.0025,col='white')

ggplot(data)+
  geom_boxplot(aes(x="DO", y=value))

ggplot(data=d2)+
  geom_violin(aes(x=Distribution,y=Value),fill='grey',trim=F)+
  geom_segment(aes(
    x=match(Distribution,levels(Distribution))-0.1,
    xend=match(Distribution,levels(Distribution))+0.1,
    y=Value,yend=Value),
               col='black'
  )

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
#                                 REQUIRES DATA FROM MARK ANALYSES!
#             #FIXME: Should be moved since it depends on other results and is thus out of order

#granivore data from dissertation chapter table 2-3
Psi = c(0.09, 0.13, 0.08, 0.14, 0.48, 0.23, 0.49, 0.41)
S = c(0.76, 0.78, 0.80, 0.81, 0.67, 0.81, 0.76, 0.62)
distbench = c(28.36, 32.64, 25.57, 36.22, 93.05, 74.29, 29.12, 93.30)
littersize = c(2.37, 2, 3.6, 4.72, 2.53, 3.6, 4, 4.29) #from Mammals of Arizona - Hoffmeister

lm1 = lm(S~Psi)
lm2 = lm(S~distbench)
lm3 = lm(littersize~S)
lm4 = lm(littersize~distbench)

# S vs. Psi
plot(Psi, S, pch = 19, xlim = c(0:1), ylim = c(0:1), bty = "n", cex.axis = 1.5, cex.lab = 1.5, cex = 1.5,
     ylab = "survival probability", xlab = "long distance movement probability")
    abline(lm1, col = "turquoise")
# S vs. distance benchmark
plot(distbench, S, pch = 19, bty = "n", xlab = "transition distance (meters)", ylab = "survival", 
     xlim = c(0,100), ylim = c(0,1), cex.axis = 1.5, cex.lab = 1.5, cex = 1.5)
    abline(lm2, col = "turquoise")
# littersize vs. S
plot(S, littersize, pch = 19, bty = "n", xlab = "survival probability", ylab = "mean litter size",
     xlim = c(0,1), ylim = c(1,6))
#littersize vs. distance benchmark
plot(distbench, littersize, pch = 19)

  