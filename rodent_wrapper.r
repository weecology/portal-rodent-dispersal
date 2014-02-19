# Code for working with individual-level rodent data
# movement with stakes

library(reshape2)
library(ggplot2)
library(calibrate)
library(fields)
library(stringr)
library(plyr)
library(gridExtra)
library(ggbiplot)
library(stringr)

#---------------------------------------------------------------------------------
#          setup - select wd, import data, source code,  file to collect results
#---------------------------------------------------------------------------------
#set working directory
wd = "/Users/sarah/Documents/GitHub/portal-rodent-dispersal"
#wd = "C:/Users/sarah/Documents/GitHub/portal-rodent-dispersal"
setwd(wd)
source("movement_fxns.R")

#import all data
all = read.table('rawdata/all_1980-2009.txt', sep = ',', header = TRUE)

#import species information
species_table = read.csv('rawdata/species.csv', sep = ',', header = TRUE)

#import cleaned data, if next step (clean up the data) has previously been run
allclean = read.csv('rawdata/cleaned_1989-2009.csv', sep = ',', header = TRUE)
  #all7 == allclean, if have allclean, can skip the datacleaning steps
  all7=allclean

#---------------------------------------------------------------------------------
#          clean up the data -- skip this section if allclean exists
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

# make sure when note2 == "*" it is counted as a new tag
# necessary if using ALL data (ear and toe tags)
# returns the dataset with new IDs for checking for duplicate tags that occur over a suspiciously long time period
tags = unique(all2$tag)
all3 = starred_tags(all2, tags, 9, 16)

#check for dead individuals, give all with same tag after marked dead a new unique ID
tags = unique(all3$tag)
all4 = is_dead(all3, tags, 9, 16)

# check for individuals with same tag, but captured over long timespan (may be able to separate them) 
# necessary if using ALL data (ear and toe tags)
# returns the dataset with new IDs for those that can easily be sorted out based on sex and year
tags = unique(all4$tag)
dups = is_duplicate_tag(all4, tags, 10, 9, 16) #check to see if can separate tags based

#eliminate bad data based on tags identified in dups$bad
duptags = unique(dups$bad$tag)
all5 = dups$data[-which(dups$data$tag %in% duptags),] #delete rows flagged as duplicates without clear resolution

tags = unique(all5$tag)
same = same_period(all5, tags)

#eliminate tags that appear more than once in the same period - questionable data
sametags = unique(same$tag)
all6 = all5[-which(all5$tag %in% sametags),]

# get rid of 'bad data'; deletes data where species is inconsistent. 
all7 = subsetDat(all6)

#subset data for years of analysis
all7 = subset(all7, yr > 1988)


#---------------------------------------------------------------------------------
#          calculate life-history details - temporal persistence, reproduction
#---------------------------------------------------------------------------------
spplist = unique(all7$species)
totalyears = c(1989:2009) #will need to be adjusted based on time frame want to use
controlplots = c(1,2,4,8,9,11,12,14,17,22)
months = count_months(all7, totalyears)

persistence = data.frame(species=NA, propyrs=NA, propmos=NA, meanabun=NA, maxabun=NA, reprod=NA, bodysize=NA)
  outcount = 1
yearly_control_abundance = data.frame(year = totalyears)
reprod_mo_yr = data.frame(year=NA, month=NA, proprepro=NA, numfemales=NA, species=NA)
avg_mo_reprod = data.frame(species=NA, month=NA, proprepro=NA)

for (i in 1:length(spplist)){
  #subset species data
  spdata = subset(all7, species == spplist[i])
  
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
  persistence[outcount,1] = as.character(spplist[i])
  persistence[outcount,2:5] = c(round(conpropyrs,4), round(conavgmos,4), round(meanconabun,4), maxconabun)
  yearly_control_abundance = cbind(yearly_control_abundance, conabun)
  
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
      persistence[outcount,6] = round(mean(irep$num_reprod),2)
  
  # record mean body size for individuals captured for the species (excluding juveniles and pregnant individuals)
  massdata = spdata[which(spdata$pregnant != "P" & spdata$reprod != "J"),]
  persistence[outcount,7] = round(mean(massdata$wgt, na.rm=TRUE),2)
  
  outcount = outcount + 1
}

#add species names to the dataframe of abundance vectors
names(yearly_control_abundance) = c("year", as.character(spplist))
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


#Identify Core species based on annual persistence, following Coyle et al. 2013 (American Naturalist)
corespecies = persistence[which(persistence$propyrs >= 0.666),1]
intermediatespecies = persistence[which(persistence$propyrs > 0.333 & persistence$propyrs < 0.666),1]
transientspecies = persistence[which(persistence$propyrs <= 0.333),1]

#Categorize species by feeding guild
granivores = c("DO", "DM", "DS", "PB", "PP", "PF", "PH", "PI",
               "PE", "PM", "PL", "RM", "RF", "RO", "BA")
folivores = c("SH", "SF", "SO", "NAO")
carnivores = c("OT", "OL")

#add columns to persistence (for later plotting) for status and guild
persistence$guild = rep(NA, nrow(persistence))
persistence$status = rep(NA, nrow(persistence))

for (row in 1:nrow(persistence)){
  if (persistence[row,]$species %in% granivores){ persistence[row,]$guild = "granivore" }
  else if (persistence[row,]$species %in% folivores) { persistence[row,]$guild = "folivore" }
  else { persistence[row,]$guild = "carnivore" }

  if (persistence[row,]$species %in% corespecies){ persistence[row,]$status = "core" }
  else if (persistence[row,]$species %in% intermediatespecies) { persistence[row,]$status = "intermediate" }
  else { persistence[row,]$status = "transient" }
}

#add status and guild to yrcontrolabuns
yrcontrolabuns$guild = rep(NA, nrow(yrcontrolabuns))
yrcontrolabuns$status = rep(NA, nrow(yrcontrolabuns))

for (row in 1:nrow(yrcontrolabuns)){
  if (yrcontrolabuns[row,]$species %in% granivores){ yrcontrolabuns[row,]$guild = "granivore" }
  else if (yrcontrolabuns[row,]$species %in% folivores) { yrcontrolabuns[row,]$guild = "folivore" }
  else { yrcontrolabuns[row,]$guild = "carnivore" }
  
  if (yrcontrolabuns[row,]$species %in% corespecies){ yrcontrolabuns[row,]$status = "core" }
  else if (yrcontrolabuns[row,]$species %in% intermediatespecies) { yrcontrolabuns[row,]$status = "intermediate" }
  else { yrcontrolabuns[row,]$status = "transient" }
}


#---------------------------------------------------------------------------------
#          calculate movement distances, multi-state capture histories
#---------------------------------------------------------------------------------

#set names for the distance list
spnames = as.character(spplist)

#create lists with sublist names reflecting each species
meterlist = sapply(spnames,function(x) NULL)
taglist = sapply(spnames,function(x) NULL)


for (i in 1:length(spplist)){
  #subset species data
  spdata = subset(all7, species == spplist[i])
  
    # get a vector unique tags, then get a vector of distances moved for all recaptured individuals, by SPECIES
    tags = unique(spdata$tag)
    print (paste(spplist[i], length(tags), sep = " "))
    print (paste(spplist[i], round(mean(spdata$wgt, na.rm=T),2), sep = " "))
    mtrs = distance_moved(spdata, tags)
    meterlist[i] = list(mtrs)
    taglist[i] = list(tags)
}

#Identify breakpoint, modal, mean and max distance for all species in meterlist
breakpoint = sapply(meterlist,function(x){
  expm1(mean(log1p(x)) + sd(log1p(x)))
})

modal_distance = sapply(meterlist,function(x){
  as.numeric(names(sort(-table(x)))[1])
})

mean_distance <-sapply(meterlist,function(x){
  as.numeric(mean(x))
})

max_distance <-sapply(meterlist,function(x){
  as.numeric(max(x))
})

persistence = cbind(persistence, breakpoint, modal_distance, mean_distance, max_distance)


#---------------------------------- From this point on, the analysis will be separated by foraging guild
#------------------------ Granivores, Folivores, and Carnivores

# concatenate core guild data - used to ask if these species behave differently from others
# Check that the definition of "core" is still the same, need to change this chunk of code by hand, if necessary
coremeters = meterlist[which(names(meterlist) %in% corespecies)]
coregran = unlist(coremeters[which(names(coremeters) %in% granivores)], use.names=F) 
corefoli = unlist(coremeters[which(names(coremeters) %in% folivores)], use.names=F)
corecarn = unlist(coremeters[which(names(coremeters) %in% carnivores)], use.names=F)

  # find breakpoints to use in MARK data structure for future analyses
  # data reasonably well fits a lognormal distribution (eyeball and J. Powell)
  # breakpoint = mean(logdata) + sd(logdata) of all the distances traveled by recaptured individuals    
  # using log1p, and back transforming using expm1 should solve the problem of lots of zeros 
  coregran_brkpt = expm1(mean(log1p(coregran)) + sd(log1p(coregran)))
  corefoli_brkpt = expm1(mean(log1p(corefoli)) + sd(log1p(corefoli)))
  corecarn_brkpt = expm1(mean(log1p(corecarn)) + sd(log1p(corecarn)))

#concatenate data for granivores only, and 
#add in the transition, modal, mean, and max distances traveled by each species for later plotting
graniv_persist = persistence[which(persistence$species %in% granivores),] 
graniv_persist$oneval = graniv_persist$propyrs * graniv_persist$propmos


#----------------------------------------------------------------------
#             Get MARK capture histories for granivores
#----------------------------------------------------------------------
spplist = granivores
periods = c(130:380) # all sampling periods 1989-2009
all_excl = c(5, 7, 10, 16, 23, 24) 
krat_excl = c(5, 7, 10, 16, 23, 24, 3, 6, 13, 15, 18, 19, 20, 21)

for (i in 1:length(spplist)){
  
  #subset species data, for each species in turn
  spdata = subset(all7, species == spplist[i])
  
  if (spplist[i] %in% list("DM", "DS", "DO")) { exclosures = krat_excl} 
  else { exclosures = all_excl} 
  
  #the first species begins the new data matrix for MARK
  if (i == 1) {
  MARK = noplacelikehome(spdata, periods, exclosures, coregran_brkpt)
  }
  
  #all subsequent species are appended onto the end of the existing MARK data matrix
  else {
  nextMARK = noplacelikehome(spdata, periods, exclosures, coregran_brkpt)
  MARK = rbind(MARK, nextMARK)
  }
}

MARK = as.data.frame(MARK[,1:3])
names(MARK) = c("ch", "freq", "species")

#separate into files based on species
domark = MARK[which(MARK[,3]=="DO"),]
dmmark = MARK[which(MARK[,3]=="DM"),]
dsmark = MARK[which(MARK[,3]=="DS"),]
pbmark = MARK[which(MARK[,3]=="PB"),]
ppmark = MARK[which(MARK[,3]=="PP"),]
pfmark = MARK[which(MARK[,3]=="PF"),]
pemark = MARK[which(MARK[,3]=="PE"),]
pmmark = MARK[which(MARK[,3]=="PM"),]
rmmark = MARK[which(MARK[,3]=="RM"),]
transmark = MARK[which(MARK[,3] %in% transientspecies),]
  transmark[,3] = "TR"
  levels(transmark[,3]) = "TR"

write.table(domark, file = "mark_datafiles//do_mark.txt", sep=" ", row.names = F)
write.table(dmmark, file = "mark_datafiles//dm_mark.txt", sep=" ", row.names = F)
write.table(dsmark, file = "mark_datafiles//ds_mark.txt", sep=" ", row.names = F)
write.table(pbmark, file = "mark_datafiles//pb_mark.txt", sep=" ", row.names = F)
write.table(ppmark, file = "mark_datafiles//pp_mark.txt", sep=" ", row.names = F)
write.table(pfmark, file = "mark_datafiles//pf_mark.txt", sep=" ", row.names = F)
write.table(pemark, file = "mark_datafiles//pe_mark.txt", sep=" ", row.names = F)
write.table(pmmark, file = "mark_datafiles//pm_mark.txt", sep=" ", row.names = F)
write.table(rmmark, file = "mark_datafiles//rm_mark.txt", sep=" ", row.names = F)
write.table(transmark, file = "mark_datafiles//trans_mark.txt", sep=" ", row.names = F)

#write.table(MARK, file = "mark_datafiles//gran_mark.inp", row.names = F, col.names = F, quote = F)


#--------------- Get MARK capture histories for folivores
#------------------------------
spplist = folivores
periods = c(130:380) #1989-2009
all_excl = c(5, 7, 10, 16, 23, 24) 

for (i in 1:length(spplist)){
  
  #subset species data, for each species in turn
  spdata = subset(all7, species == spplist[i])
   exclosures = all_excl
  
  #the first species begins the new data matrix for MARK
  if (i == 1) {
    MARK = noplacelikehome(spdata, periods, exclosures, corefoli_brkpt)
  }
  
  #all subsequent species are appended onto the end of the existing MARK data matrix
  else {
    nextMARK = noplacelikehome(spdata, periods, exclosures, corefoli_brkpt)
    MARK = rbind(MARK, nextMARK)
  }
}

MARK = as.data.frame(MARK[,1:3])
names(MARK) = c("ch", "freq", "species")

#separate into files based on species
naomark = MARK[which(MARK[,3]=="NAO"),]
shmark = MARK[which(MARK[,3]=="SH"),]
sfmark = MARK[which(MARK[,3]=="SF"),]
somark = MARK[which(MARK[,3]=="SO"),]

write.table(naomark, file = "mark_datafiles//nao_mark.txt", sep=" ", row.names = F)
write.table(shmark, file = "mark_datafiles//sh_mark.txt", sep=" ", row.names = F)
write.table(sfmark, file = "mark_datafiles//sf_mark.txt", sep=" ", row.names = F)
write.table(somark, file = "mark_datafiles//so_mark.txt", sep=" ", row.names = F)

#write.table(MARK, file = "mark_datafiles//foli_mark.inp", row.names = F, col.names = F, quote = F)


#-----------------Get MARK capture histories for carnivores
#------------------------------
spplist = carnivores
periods = c(130:380) #1989-2009
all_excl = c(5, 7, 10, 16, 23, 24) 

for (i in 1:length(spplist)){
  
  #subset species data, for each species in turn
  spdata = subset(all7, species == spplist[i])
  exclosures = all_excl
  
  #the first species begins the new data matrix for MARK
  if (i == 1) {
    MARK = noplacelikehome(spdata, periods, exclosures, corecarn_brkpt)
  }
  
  #all subsequent species are appended onto the end of the existing MARK data matrix
  else {
    nextMARK = noplacelikehome(spdata, periods, exclosures, corecarn_brkpt)
    MARK = rbind(MARK, nextMARK)
  }
}

MARK = as.data.frame(MARK[,1:3])
names(MARK) = c("ch", "freq", "species")

#separate into files based on species
otmark = MARK[which(MARK[,3]=="OT"),]
olmark = MARK[which(MARK[,3]=="OL"),]

write.table(otmark, file = "mark_datafiles//ot_mark.txt", sep=" ", row.names = F)
write.table(olmark, file = "mark_datafiles//ol_mark.txt", sep=" ", row.names = F)

#write.table(MARK, file = "mark_datafiles//carn_mark.inp", row.names = F, col.names = F, quote = F)


#-----------------------------------------------------------------------------------
#        Make a table of total capture/recapture and gaps in data for each species
#-----------------------------------------------------------------------------------
# this code only uses the most recently generated MARK data table, 
# so would need to run separately for granivores, folivores, and carnivores 
spplist = unique(MARK[,2])
TableDat = data.frame("species"=1, "numindiv"=1, "numindivrecap"=1, "nocaps"=1, "norecaps"=1)

for (i in 1:length(spplist)){
  
  #subset species data, for each species in turn
  spdata = MARK[which(MARK[,2]==spplist[i]),]
  
  sp = spplist[i]
  numindiv = nrow(spdata)
  
  #break MARK data out into a list
  list.all = list()
  for (row in 1:nrow(spdata)){
    dat = spdata[row,1]
    list.all[row] = dat
  }
  
  #make a new matrix with each period as a value
  out<-t(sapply(list.all,function(x){
    as.numeric(as.vector(strsplit(x,"")[[1]]))
  }))
  
  #For each row find the sum things that are not 0
  recaps<-apply(out,1,function(x){
    length(which(!x %in% 0))-1
  })
  
  numinds_recaps = length(recaps[recaps>0])
  
  #count the num of captures per period
  caps = apply(out,2,function(x) {
    length(which(!x %in% 0))
  })

  #count the num of recaptures per period
  recaps2 = apply(out,2,function(x) {
    length(which(x %in% 2))
  })
  
  nummos_nocaps = length(caps[caps == 0])
  nummos_norecaps = length(recaps2[recaps2 == 0])
  
  results = c(sp, numindiv, numinds_recaps, nummos_nocaps, nummos_norecaps)
  TableDat[i,] <- results
}

print(TableDat)


#-----------------------------------------------------------------------------------
#        If MARK_analyses.r has already been run, and results saved as .csv files
#          analyze the data from the Program Mark analysis
#-----------------------------------------------------------------------------------
#---------- concatenate results
#grab all the .csv files to loop over and summarize results
rfiles = list.files(paste(getwd(), "/mark_output", sep=""), pattern = "real", full.name=T)

# loop thru the files to make a new dataframe with the estimated parameters
estimates = data.frame(species=NA, S=1, S_se=1, p=1, p_se=1, Psi=1, Psi_se = 1)
outcount = 1

for (f in 1:length(rfiles)){
  dat = read.csv(rfiles[f], header=T, sep=",")
  spname = str_sub(rfiles[f],-6,-5)
  estimates[outcount,1] = spname
  estimates[outcount,2:7] = c(dat[1,2], dat[1,3], dat[2,2], dat[2,3], dat[3,2], dat[3,3])
  outcount = outcount + 1
}

#assign all transient species the same MARK estimates (for now)
estimates2 = estimates[16,]
transgran = granivores[which(granivores %in% transientspecies)]
for (i in 1:length(transgran)){
  estimates2$species = transgran[i]
  estimates=rbind(estimates, estimates2)
}
#change "AO" to "NAO" to match naming schema - shortened in MARK_analyses.r
estimates[which(estimates$species == "AO"),1] = "NAO" 

# merge GRANIVORE data for analysis and pca plots 
mgran = merge(graniv_persist, estimates)

pcagraniv = mgran[,c(1, 2, 5, 6, 7, 10, 11, 15, 17, 19)]
catgraniv = mgran[,c(1, 8, 9)]
rownames(pcagraniv) = pcagraniv$species
pcagraniv = pcagraniv[,-1]

# merge ALL data for analysis and pca plots
mall = merge(persistence, estimates, by = c("species", "species"))
mall = merge(mall, species_table)

pcadat = mall[,c(1, 2, 5, 6, 7, 10, 11, 14, 16, 18)]
catdat = mall[,c(1, 8, 9, 20)]
rownames(pcadat) = pcadat$species
pcadat = pcadat[,-1]


#---------------------------------------------------------------------------------
#                                  plot results
#---------------------------------------------------------------------------------

#------------------------ PCA biplot of species traits and estimates

# Standardize the matrix to correct for different units by subtracting the
# means and dividing by sd
zscore <- apply(pcagraniv, 2, function(x) {
  y <- (x - mean(x))/sd(x)
  return(y)
  })
rownames(zscore) <- rownames(pcagraniv)

#run pca analysis
trait_pc = prcomp(zscore)

#Use ggplot biplot to color by grouped categories
toCol = catgraniv[catgraniv$species %in% rownames(trait_pc$x),"status"]

#Label species names and clades, circles cover normal distribuiton of groups
ggbiplot(trait_pc, groups=toCol, labels=rownames(trait_pc$x), label.size = 3, varname.size = 4) + theme_classic() +
  theme(text = element_text(size=20))


#------------------------ PCA biplot for all the species

# Standardize the matrix to correct for different units by subtracting the
# means and dividing by sd
zscore <- apply(pcadat, 2, function(x) {
  y <- (x - mean(x))/sd(x)
  return(y)
  })
rownames(zscore) <- rownames(pcadat)

#run pca analysis
trait_pc<-prcomp(zscore)

#Use dev libary to ggplot PCA, color by clades
#Try the ggplot biplot to color by clades (or later, behavioral roles)
toCol = catdat[catdat$species %in% rownames(trait_pc$x),"status"]
toColGuild = catdat[catdat$species %in% rownames(trait_pc$x),"guild"]
toColFam = catdat[catdat$species %in% rownames(trait_pc$x),"family"]

#Label species names and clades, ellipses cover normal distribuiton of temporal groups
ggbiplot(trait_pc, groups=toCol, labels=rownames(trait_pc$x), ellipse=TRUE, label.size = 3, varname.size = 4) + theme_classic() +
  theme(text = element_text(size=20))

#Label species names and clades, circles cover normal distribuiton of guilds
ggbiplot(trait_pc, groups=toColGuild, labels=rownames(trait_pc$x), ellipse=TRUE, label.size = 3, varname.size = 4) + theme_classic() +
  theme(text = element_text(size=20))

#Label species names and clades, circles cover normal distribuiton of families
ggbiplot(trait_pc, groups=toColFam, labels=rownames(trait_pc$x), ellipse=TRUE, label.size = 3, varname.size = 4) + theme_classic() +
  theme(text = element_text(size=20))


#------------------------ plot the Mark results for all the species

# Use only the real estimates, where transient granivores were lumped
est = estimates[c(1:16),]

# Add status 
toStatus = catdat[catdat$species %in% est$species,"status"]
toStatus = append(toStatus, "transient") #since TR is a group instead of an actual "species"

est = cbind(est,toStatus)

SbyPsi = ggplot(est, aes(Psi, S, col=toStatus)) + geom_point(size = 3) + theme_classic() + 
  theme(text = element_text(size=20)) + scale_colour_hue(guide = "none") +
  xlab("Long-distance movement probability") + ylab("Survival probability") + 
  geom_errorbar(aes(x = Psi, ymin = S - S_se, ymax = S + S_se), width=0.01) +
  geom_errorbarh(aes(xmin = Psi - Psi_se, xmax = Psi + Psi_se))

Sbyp = ggplot(est, aes(p, S, col=toStatus)) + geom_point(size = 3) + theme_classic() + 
  theme(text = element_text(size=20)) + scale_colour_hue(guide = "none")+
  xlab("recapture probability") + ylab("Survival probability") + 
  geom_errorbar(aes(x = p, ymin = S - S_se, ymax = S + S_se), width=0.01) +
  geom_errorbarh(aes(xmin = p - p_se, xmax = p + p_se)) + ggtitle("all species")

Psibyp =  ggplot(est, aes(Psi, p, col=toStatus)) + geom_point(size = 3) + theme_classic() + 
  theme(text = element_text(size=20)) + scale_colour_hue(guide = "none") +
  xlab("long-distance movement probability") + ylab("recapture probability") + 
  geom_errorbar(aes(x = Psi, ymin = p - p_se, ymax = p + p_se), width=0.01) +
  geom_errorbarh(aes(xmin = Psi - Psi_se, xmax = Psi + Psi_se))

grid.arrange(SbyPsi, Sbyp, Psibyp, nrow=1)


#------------------------ plot abundance vs. years, ala core v. transient literature

status_plot = ggplot(persistence, aes(propyrs, propmos)) + geom_point(aes(shape=guild), size = 3) + 
  theme_classic() + xlab("proportion years present") +
  ylab("proportion of months present") + xlim(0,1) + ylim(0,1) + 
  geom_vline(xintercept=0.66, linetype="dotted", col = "red", size = 1.5) + 
  geom_vline(xintercept=0.33, linetype="dotted", col = "red", size = 1.5) + 
  ggtitle("Rodents 1989 - 2009") + geom_text(aes(label = species), hjust=0, vjust=0) +
  theme(text = element_text(size=20))


#------------------------- plot abundance for all species across timeseries
ggplot(yrcontrolabuns, aes(x=year, y=abun, group=species)) + 
  geom_line(aes(col=status), size=1.5) + theme_bw() +
  theme(text = element_text(size=20))


#------------------------- plot monthly reproduction
ggplot(avg_mo_reprod, aes(month, proprepro)) + geom_point() + theme_bw() +
  geom_line() + facet_wrap(~species)


#------------------------- plot meters traveled by all species
#plot modal distance by persistence for all granivores, color code points by status
modal_dist = ggplot(graniv_persist, aes(propyrs, modal_distance)) + theme_classic() +
  geom_point(aes(col=as.factor(status), shape = as.factor(status)), size=5)  + 
  xlab("proportion of years present") + ylab("Modal Distance between trap locations") +
  geom_text(aes(label=species), hjust=0, vjust=0) +
  theme(text = element_text(size=20))

#plot modal distance by persistence for all granivores, color code points by status
modal_all = ggplot(persistence, aes(propyrs, modal_distance)) + theme_classic() +
  geom_point(aes(col=as.factor(status), shape = as.factor(status)), size=5)  + 
  xlab("proportion of years present") + ylab("Modal Distance between trap locations") +
#  geom_text(aes(label=species), hjust=0, vjust=0) +
  theme(text = element_text(size=20))

grid.arrange(modal_dist, modal_all, nrow = 1)


#------------------------ plot histograms of each of the groups (by status) movements

core = data.frame(core = c(meterlist$DO, meterlist$DM, meterlist$RM, meterlist$PE, meterlist$PP, meterlist$PB, meterlist$PF))
interm = data.frame(interm = c(meterlist$PM, meterlist$DS))
trans = data.frame(trans=c(meterlist$RF, meterlist$BA, meterlist$PH, meterlist$PI, meterlist$PL, meterlist$RO))

COREplot = ggplot(core, aes(core)) + geom_histogram() + theme_bw() + 
  theme(text = element_text(size=20)) + xlab("distance (m) between recaptures") + ggtitle("Core") +
  xlim(0,550)

INTERplot = ggplot(interm, aes(interm)) + geom_histogram() + theme_bw() + 
  theme(text = element_text(size=20)) + xlab("distance (m) between recaptures") + ggtitle("Intermediate") +
  xlim(0,550)

TRANSplot = ggplot(trans, aes(trans)) + geom_histogram() + theme_bw() + 
  theme(text = element_text(size=20)) + xlab("distance (m) between recaptures") + ggtitle("Transient") +
  xlim(0,550)

grid.arrange(COREplot, INTERplot, TRANSplot, nrow=1)

#------------------------ plot histograms of each of the species movements
splist = names(meterlist)
breaks = c(0, seq(3,600, by=6))
i = 1
plot = list() 

for (i in 1:length(splist)){
  df = data.frame(meterlist[i])
  names(df) = "spp"
  plot[[i]] = ggplot(df, aes(spp)) + geom_histogram(breaks=breaks, position="dodge") + theme_classic() + ggtitle(splist[i]) +
    theme(text = element_text(size=20)) + xlab("recapture distance")
  i = i + 1
}

#core histograms
grid.arrange(plot[[3]], plot[[10]], plot[[21]], plot[[12]], plot[[13]], 
             plot[[8]], plot[[4]], plot[[7]], plot[[1]], plot[[11]], ncol = 4)
#intermediate histograms
grid.arrange(plot[[9]], plot[[2]], plot[[5]], plot[[16]], ncol = 4)
#transient histograms
grid.arrange(plot[[15]], plot[[19]], plot[[20]], plot[[14]], plot[[6]], 
             plot[[17]], plot[[17]], ncol = 4)


#-------------------------- Comparing the CMR analysis estimates
lm1 = lm(S ~ Psi, data = mall)
PsiS = ggplot(mall, aes(Psi, S)) + geom_point(size = 2) + stat_smooth(method = "lm") + theme_classic() + 
  theme(text = element_text(size=20)) + ggtitle(paste("r2 =", round(summary(lm1)$r.squared,2)))

lm1 = lm(p ~ Psi, data = mall)
Psip = ggplot(mall, aes(Psi, p)) + geom_point(size = 2) + stat_smooth(method = "lm") + theme_classic() + 
  theme(text = element_text(size=20)) + ggtitle(paste("r2 =", round(summary(lm1)$r.squared,2)))

lm1 = lm(Psi ~ propyrs, data = mall)
yrsPsi = ggplot(mall, aes(propyrs, Psi)) + geom_point(size = 2) + stat_smooth(method = "lm") + theme_classic() + 
  theme(text = element_text(size=20)) + ggtitle(paste("r2 =", round(summary(lm1)$r.squared,2)))

lm1 = lm(S ~ propyrs, data = mall)
yrsS = ggplot(mall, aes(propyrs, S)) + geom_point(size = 2) + stat_smooth(method = "lm") + theme_classic() + 
  theme(text = element_text(size=20)) + ggtitle(paste("r2 =", round(summary(lm1)$r.squared,2)))

lm1 = lm(modal_distance ~ propyrs, data = mall)
yrsdist = ggplot(mall, aes(propyrs, modal_distance)) + geom_point(size = 2) + stat_smooth(method = "lm") + theme_classic() + 
  theme(text = element_text(size=20)) + ggtitle(paste("r2 =", round(summary(lm1)$r.squared,2)))

lm1 = lm(maxabun ~ propyrs, data = mall)
yrsabun = ggplot(mall, aes(propyrs, maxabun)) + geom_point(size = 2) + stat_smooth(method = "lm") + theme_classic() + 
  theme(text = element_text(size=20)) + ggtitle(paste("r2 =", round(summary(lm1)$r.squared,2)))

lm1 = lm(meanabun ~ propyrs, data = mall)
yrsmean = ggplot(mall, aes(propyrs, meanabun)) + geom_point(size = 2) + stat_smooth(method = "lm") + theme_classic() + 
  theme(text = element_text(size=20)) + ggtitle(paste("r2 =", round(summary(lm1)$r.squared,2)))

yrsreprod = ggplot(mall, aes(propyrs, reprod)) + geom_point(size = 2) + theme_classic() + 
  theme(text = element_text(size=20))

yrssd = ggplot(mall, aes(propyrs, breakpoint)) + geom_point(size = 2) + theme_classic() + 
  theme(text = element_text(size=20))

grid.arrange(PsiS, Psip, yrsPsi, yrsS, yrsdist, yrssd, yrsabun, yrsmean, yrsreprod, nrow = 3)


# other plots
ggplot(mall, aes(Psi, reprod)) + geom_point(size = 2) + stat_smooth(method = "lm") + theme_classic() + 
  theme(text = element_text(size=20))

ggplot(mall, aes(bodysize, Psi)) + geom_point(size = 2) + stat_smooth(method = "lm") + theme_classic() + 
  theme(text = element_text(size=20))

ggplot(mall, aes(bodysize, S)) + geom_point(size = 2) + stat_smooth(method = "lm") + theme_classic() + 
  theme(text = element_text(size=20))


#-------------------------- Output summary info for all species
core_m = mall[which(mall$species %in% corespecies),]
inter_m = mall[which(mall$species %in% intermediatespecies),]
trans_m = mall[which(mall$species %in% transientspecies),]

mean(core_m$breakpoint)
mean(inter_m$breakpoint)
mean(trans_m$breakpoint)

mean(core_m$Psi)
mean(inter_m$Psi)
mean(trans_m$Psi)

mean(core_m$p)
mean(inter_m$p)
mean(trans_m$p)

mean(core_m$S)
mean(inter_m$S)
mean(trans_m$S)

#-------------------------- Output summary info for GRANIVORES only
core_m = mgran[which(mgran$species %in% corespecies),]
inter_m = mgran[which(mgran$species %in% intermediatespecies),]
trans_m = mgran[which(mgran$species %in% transientspecies),]

mean(core_m$breakpoint)
mean(inter_m$breakpoint)
mean(trans_m$breakpoint)

mean(core_m$Psi)
mean(inter_m$Psi)
mean(trans_m$Psi)

mean(core_m$p)
mean(inter_m$p)
mean(trans_m$p)

mean(core_m$S)
mean(inter_m$S)
mean(trans_m$S)


#-------------------------- Reproduction Figures
pdf("reproductive_events_per_year.pdf")
par(mfrow=c(3,3))

#plot the relative abundance of females who represent each number of reproductive events per year
plot(NA, NA, xlim = c(0, 4), xaxp = c(0, 4, 4), xlab = "number reproductive events per year", 
     ylab = "proportion females", type = "b", pch = 19, ylim = c(0, 1), yaxp = c(0, 1, 4), bty = "n", 
     main = "core granivores", cex.axis = 1.5, cex.lab = 1.5)
points(c(0:3), table(DOirep$num_reprod)/sum(table(DOirep$num_reprod)), type = "b", pch = 19, col = "black")
points(c(0:4), table(DMirep$num_reprod)/sum(table(DMirep$num_reprod)), type = "b", pch = 19, col = "black")
points(c(0:4), table(PBirep$num_reprod)/sum(table(PBirep$num_reprod)), type = "b", pch = 19, col = "black")
points(c(0:4), table(PPirep$num_reprod)/sum(table(PPirep$num_reprod)), type = "b", pch = 19, col = "black")
points(c(0:3), table(RMirep$num_reprod)/sum(table(RMirep$num_reprod)), type = "b", pch = 19, col = "black")
points(c(0:2), table(PEirep$num_reprod)/sum(table(PEirep$num_reprod)), type = "b", pch = 19, col = "black")
points(c(0:4), table(PFirep$num_reprod)/sum(table(PFirep$num_reprod)), type = "b", pch = 19, col = "black")

plot(NA, NA, xlim = c(0, 4), xaxp = c(0, 4, 4), xlab = "number reproductive events per year", 
     ylab = "proportion females", type = "b", pch = 19, ylim = c(0, 1), yaxp = c(0, 1, 4), bty = "n", 
     main = "intermediate granivores", cex.axis = 1.5, cex.lab = 1.5)
points(c(0:2), table(PMirep$num_reprod)/sum(table(PMirep$num_reprod)), type = "b", pch = 15, col = "black")
points(c(0,1,3), table(DSirep$num_reprod)/sum(table(DSirep$num_reprod)), type = "b", pch = 15, col = "black")

plot(NA, NA, xlim = c(0, 4), xaxp = c(0, 4, 4), xlab = "number reproductive events per year", 
     ylab = "proportion females", type = "b", pch = 19, ylim = c(0, 1), yaxp = c(0, 1, 4), bty = "n", 
     main = "transient granivores", cex.axis = 1.5, cex.lab = 1.5)
points(c(0:1), table(PLirep$num_reprod)/sum(table(PLirep$num_reprod)), type = "b", pch = 15, col = "black")
points(c(0:1), table(ROirep$num_reprod)/sum(table(ROirep$num_reprod)), type = "b", pch = 15, col = "black")
points(0, table(RFirep$num_reprod)/sum(table(RFirep$num_reprod)), type = "b", pch = 15, col = "black")
points(c(0:1), table(PIirep$num_reprod)/sum(table(PIirep$num_reprod)), type = "b", pch = 15, col = "black")
points(c(0:1), table(PHirep$num_reprod)/sum(table(PHirep$num_reprod)), type = "b", pch = 15, col = "black")
points(c(0:2), table(BAirep$num_reprod)/sum(table(BAirep$num_reprod)), type = "b", pch = 15, col = "black")

# For Folivores:
plot(NA, NA, xlim = c(0, 4), xaxp = c(0, 4, 4), xlab = "number reproductive events per year", 
     ylab = "proportion females", type = "b", pch = 19, ylim = c(0, 1), yaxp = c(0, 1, 4), bty = "n", 
     main = "core folivores", cex.axis = 1.5, cex.lab = 1.5)
points(c(0:3), table(NAOirep$num_reprod)/sum(table(NAOirep$num_reprod)), type = "b", pch = 15, col = "black")

plot(NA, NA, xlim = c(0, 4), xaxp = c(0, 4, 4), xlab = "number reproductive events per year", 
     ylab = "proportion females", type = "b", pch = 19, ylim = c(0, 1), yaxp = c(0, 1, 4), bty = "n", 
     main = "intermediate folivores", cex.axis = 1.5, cex.lab = 1.5)
points(c(0:1), table(SHirep$num_reprod)/sum(table(SHirep$num_reprod)), type = "b", pch = 15, col = "black")
points(c(0:1), table(SFirep$num_reprod)/sum(table(SFirep$num_reprod)), type = "b", pch = 15, col = "black")

plot(NA, NA, xlim = c(0, 4), xaxp = c(0, 4, 4), xlab = "number reproductive events per year", 
     ylab = "proportion females", type = "b", pch = 19, ylim = c(0, 1), yaxp = c(0, 1, 4), bty = "n", 
     main = "transient folivores", cex.axis = 1.5, cex.lab = 1.5)
points(c(0:2), table(SOirep$num_reprod)/sum(table(SOirep$num_reprod)), type = "b", pch = 15, col = "black")

#For Carnivores
plot(NA, NA, xlim = c(0, 4), xaxp = c(0, 4, 4), xlab = "number reproductive events per year", 
     ylab = "proportion females", type = "b", pch = 19, ylim = c(0, 1), yaxp = c(0, 1, 4), bty = "n", 
     main = "core carnivores", cex.axis = 1.5, cex.lab = 1.5)
points(c(0:3), table(OTirep$num_reprod)/sum(table(OTirep$num_reprod)), type = "b", pch = 15, col = "black")
points(c(0:2), table(OLirep$num_reprod)/sum(table(OLirep$num_reprod)), type = "b", pch = 15, col = "black")

dev.off()


# #------------------------------------------ FIGURE - for ESA talk 
# #                                 REQUIRES DATA FROM MARK ANALYSES!
# 
# #granivore data from dissertation chapter table 2-3 #FIXME - should use data from estimates (above)
# Psi = c(0.09, 0.13, 0.08, 0.14, 0.48, 0.23, 0.49, 0.41)
# S = c(0.76, 0.78, 0.80, 0.81, 0.67, 0.81, 0.76, 0.62)
# distbench = c(28.36, 32.64, 25.57, 36.22, 93.05, 74.29, 29.12, 93.30)
# littersize = c(2.37, 2, 3.6, 4.72, 2.53, 3.6, 4, 4.29) #from Mammals of Arizona - Hoffmeister
# 
# lm1 = lm(S~Psi)
# lm2 = lm(S~distbench)
# lm3 = lm(littersize~S)
# lm4 = lm(littersize~distbench)
# 
# # S vs. Psi
# plot(Psi, S, pch = 19, xlim = c(0:1), ylim = c(0:1), bty = "n", cex.axis = 1.5, cex.lab = 1.5, cex = 1.5,
#      ylab = "survival probability", xlab = "long distance movement probability")
#     abline(lm1, col = "turquoise")
# # S vs. distance benchmark
# plot(distbench, S, pch = 19, bty = "n", xlab = "transition distance (meters)", ylab = "survival", 
#      xlim = c(0,100), ylim = c(0,1), cex.axis = 1.5, cex.lab = 1.5, cex = 1.5)
#     abline(lm2, col = "turquoise")
# # littersize vs. S
# plot(S, littersize, pch = 19, bty = "n", xlab = "survival probability", ylab = "mean litter size",
#      xlim = c(0,1), ylim = c(1,6))
# #littersize vs. distance benchmark
# plot(distbench, littersize, pch = 19)
# 
#   