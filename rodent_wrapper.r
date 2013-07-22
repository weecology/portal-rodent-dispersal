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
      assign(paste0(spplist[i], 'yr'), propyrs)
    # average number of months they were seen in during years in which they were present
    avgmos = mean_win_yr_occ(spdata)
      assign(paste0(spplist[i], 'avgmos'), avgmos)
    # mean abundance within all years 
    avgabun = allyrs_abun(spdata)
      assign(paste0(spplist[i], 'avgabun'), avgabun)
    #mean abundances within all years on control plots only
    conabun = allyrs_abun(subset(spdata, plot %in% c(1,2,4,8,9,11,12,14,17,22)))
      assign(paste0(spplist[i], 'conabun'), conabun)
  
  #subset females for each species
  spdataF = subset(spdata, sex == "F")
  
    #average proportion of reproductive females by month across all years
    reprd = mean_mo_repro(spdataF) #vector with 12 items
      assign(paste0(spplist[i], 'reprd'), reprd)
    # proportion of reproductive females by month and year
    reprdyr = mo_repro(spdataF) #matrix with 120 x 5
      assign(paste0(spplist[i], 'reprdyr'), reprdyr)
    # track the number of times individual females uniquely reproduce within years
    irep = indiv_repro(spdataF) #matix with cols "tag", "year" and "num_reprod"
      assign(paste0(spplist[i], 'irep'), irep)
}

abuns = cbind(*abun)
#---------------------------------------------------------------------------------
#          calculate movement distances, multi-state capture histories
#---------------------------------------------------------------------------------

for (i in 1:length(spplist)){
  #subset species data
  spdata = subset(all2, species == spplist[i])
  
    # get a vector unique tags, then get a vector of distances moved for all recaptured individuals, by SPECIES
    tags = unique(spdata$tag)
      assign(paste0(spplist[i], 'propyrs'), tags)
    meters = distance_moved(spdata, tags)
     assign(paste0(spplist[i], 'propyrs'), meters)
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


#------------------------------------------ FIGURE 2
pdf("Fig2_temporal_occ_spp_rank.pdf", 5, 5, pointsize = 10)
par(mfrow=c(1,1))

#average rank of each species 2000-2009 on control plots
#ordered: DM, DO, PB, PP, PF, PE, PM, RM, SH, SF, NA, OT, OL
conranks = c(2.9, 4, 2.6, 1, 7.875, 7.1, 9.8, 9.333, 9, 8.6, 8.5, 4.8, 10.5)

plot(conranks[5], PFyr, xlim = c(0,1), ylim = c(1,13), ylab = "species average ranked abundance", 
     xlab = "proportion of years present", pch = 19, col = "black", cex = 1.5, cex.axis = 1.5, cex.lab = 1.5)
textxy(PFyr, conranks[5], "PF", cx = 1)
points(PPyr, conranks[4], pch = 19, col = "black", cex = 1.5)
textxy(PPyr, conranks[4],"CP", cx = 1)
points(DOyr, conranks[2], pch = 19, col = "black", cex = 1.5)
textxy(DOyr, conranks[2], "DO", cx = 1)
points(DMyr, conranks[1], pch = 19, col = "black", cex = 1.5)
textxy(DMyr, conranks[1], "DM", cx = 1)
points(PByr, conranks[3], pch = 19, col = "black", cex = 1.5)
textxy(PByr, conranks[3], "CB", cx = 1)
points(PEyr, conranks[6], pch = 19, col = "gray30", cex = 1.5)
textxy(PEyr, conranks[6], "PE", cx = 1)
points(PMyr, conranks[7], pch = 19, col = "gray30", cex = 1.5)
textxy(PMyr, conranks[7], "PM", cx = 1)
points(RMyr, conranks[8], pch = 19, col = "gray30", cex = 1.5)
textxy(RMyr, conranks[8], "RM", cx = 1)
points(OTyr, conranks[12], pch = 6, col = "gray30", cex = 1.5)
textxy(OTyr, conranks[12], "OT", cx = 1)
points(OLyr, conranks[13], pch = 6, col = "gray30", cex = 1.5)
textxy(OLyr, conranks[13], "OL", cx = 1)
points(SHyr, conranks[9], pch = 0, col = "gray30", cex = 1.5)
textxy(SHyr, conranks[9], "SH", cx = 1)
points(SFyr, conranks[10], pch = 0, col = "gray30", cex = 1.5)
textxy(SFyr, conranks[10], "SF", cx = 1)
points(NAOyr, conranks[11], pch = 0, col = "gray30", cex = 1.5)
textxy(NAOyr, conranks[11], "NA", cx = 1)

abline(v = 0.95, lty = 2, col = 'gray40', lwd = 1)
abline(h = 5, lty = 2, col = 'gray40', lwd = 1)

dev.off()

#------------------------------------------ FIGURE 3
pdf("Fig3_indiv_repro_trends.pdf", 6, 6, pointsize = 10)
par(mfrow=c(2,2))

#plot the relative abundance of females who represent each number of reproductive events per year
plot(NA, NA, xlim = c(0, 4), xaxp = c(0, 4, 4), xlab = "number reproductive events per year", 
     ylab = "proportion females", type = "b", pch = 19, ylim = c(0, 1), yaxp = c(0, 1, 4), bty = "n", 
     main = "heteromyid granivores", cex.axis = 1.5, cex.lab = 1.5)
points(c(0:3), table(DOirep$num_reprod)/sum(table(DOirep$num_reprod)), type = "b", pch = 19, col = "black")
points(c(0:4), table(DMirep$num_reprod)/sum(table(DMirep$num_reprod)), type = "b", pch = 15, col = "black")
points(c(0:4), table(PBirep$num_reprod)/sum(table(PBirep$num_reprod)), type = "b", pch = 17, col = "black")
points(c(0:4), table(PPirep$num_reprod)/sum(table(PPirep$num_reprod)), type = "b", pch = 25, col = "black")
points(c(0:2), table(PFirep$num_reprod)/sum(table(PFirep$num_reprod)), type = "b", pch = 24, col = "black")
legend('topright', c('DO', 'DM', 'CB', 'CP', 'PF'), bty = 'n', 
       col = 'black', pch = c(19,15,17,25,24), cex = 0.75)

plot(NA, NA, xlim = c(0, 4), xaxp = c(0, 4, 4), xlab = "number reproductive events per year", 
     ylab = "proportion females", type = "b", pch = 19, ylim = c(0, 1), yaxp = c(0, 1, 4), bty = "n", 
     main = "cricetid granivores", cex.axis = 1.5, cex.lab = 1.5)
points(c(0:2), table(PEirep$num_reprod)/sum(table(PEirep$num_reprod)), type = "b", pch = 19, col = "black")
points(c(0:2), table(PMirep$num_reprod)/sum(table(PMirep$num_reprod)), type = "b", pch = 15, col = "black")
points(c(0:1), table(RMirep$num_reprod)/sum(table(RMirep$num_reprod)), type = "b", pch = 17, col = "black")
legend('topright', c('PE', 'PM', 'RM'), bty = 'n', col = 'black', pch = c(19,15,17), cex = 0.75)

plot(NA, NA, xlim = c(0, 4), xaxp = c(0, 4, 4), xlab = "number reproductive events per year", 
     ylab = "proportion females", type = "b", pch = 19, ylim = c(0, 1), yaxp = c(0, 1, 4), bty = "n",
     main = "folivores", cex.axis = 1.5, cex.lab = 1.5)
points(c(0:1), table(SHirep$num_reprod)/sum(table(SHirep$num_reprod)), type = "b", pch = 19, col = "black")
points(c(0:1), table(SFirep$num_reprod)/sum(table(SFirep$num_reprod)), type = "b", pch = 15, col = "black")
points(c(0:3), table(NAOirep$num_reprod)/sum(table(NAOirep$num_reprod)), type = "b", pch = 17, col = "black")
legend('topright', c('SH', 'SF', 'NA'), bty = 'n', col = 'black', pch = c(19,15,17), cex = 0.75)

plot(NA, NA, xlim = c(0, 4), xaxp = c(0, 4, 4), xlab = "number reproductive events per year", bty = "n", 
     ylab = "proportion females", type = "b", pch = 19, ylim = c(0, 1), yaxp = c(0, 1, 4), 
     main = "carnivores", cex.axis = 1.5, cex.lab = 1.5)
points(c(0:3), table(OTirep$num_reprod)/sum(table(OTirep$num_reprod)), type = "b", pch = 19, col = "black")
points(c(0:1), table(OLirep$num_reprod)/sum(table(OLirep$num_reprod)), type = "b", pch = 15, col = "black")
legend('topright', c('OT', 'OL'), bty = 'n', col = 'black', pch = c(19,15), cex = 0.75)

dev.off()

  