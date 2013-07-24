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
all = read.table('rawdata/all_1980-2009.txt', sep = ',')

#import the data by guild
allrats = read.csv("rawdata//allrodents_1978-2012.csv") #all rodents, for temporal figure
het = read.csv("rawdata//heteromyids_2000-2009.csv")   # DO, DM, PB, PP, PF
cricet = read.csv("rawdata//cricetids_2000-2009.csv")  # PE, PM, RM
foliv = read.csv("rawdata//folivores_2000-2009.csv")   # SH, SF, NA (as NAO)
insec = read.csv("rawdata//onychomys_2000-2009.csv")   # OT, OL

# output directed to rodent_results.txt in wd. output is appended
# to existing file. output also sent to terminal. 
#sink("rodent_results.txt", append=TRUE, split=TRUE)

#---------------------------------------------------------------------------------
#          clean up the data
#---------------------------------------------------------------------------------

# change some cols from factor to character class
het$tag = as.character(het$tag); cricet$tag = as.character(cricet$tag); foliv$tag = as.character(foliv$tag); insec$tag = as.character(insec$tag)
het$species = as.character(het$species); cricet$species = as.character(cricet$species); foliv$species = as.character(foliv$species); insec$species = as.character(insec$species)
het$sex = as.character(het$sex); cricet$sex = as.character(cricet$sex); foliv$sex = as.character(foliv$sex); insec$sex = as.character(insec$sex)

# give untagged individuals a unique 7-number code
het = id_unknowns(het, 16); cricet = id_unknowns(cricet, 16); foliv = id_unknowns(foliv, 16); insec = id_unknowns(insec, 16)

# get rid of 'bad data'; deletes data that is not a pit tag, where sex is inconsistent or where species is inconsistent. 
het = subsetDat(het); cricet = subsetDat(cricet); foliv = subsetDat(foliv); insec = subsetDat(insec)

#---------------------------------------------------------------------------------
#          calculate life-history details - reproduction, temporal persistence
#---------------------------------------------------------------------------------

#average proportion of reproductive females by month across all years
doreprd = mean_mo_repro(subset(het, species == "DO" & sex == "F")); dmreprd = mean_mo_repro(subset(het, species == "DM" & sex == "F")); pfreprd = mean_mo_repro(subset(het, species == "PF" & sex == "F")); ppreprd = mean_mo_repro(subset(het, species == "PP" & sex == "F")); pbreprd = mean_mo_repro(subset(het, species == "PB" & sex == "F"));
pereprd = mean_mo_repro(subset(cricet, species == "PE" & sex == "F")); pmreprd = mean_mo_repro(subset(cricet, species == "PM" & sex == "F")); rmreprd = mean_mo_repro(subset(cricet, species == "RM" & sex == "F"))
shreprd = mean_mo_repro(subset(foliv, species == "SH" & sex == "F")); sfreprd = mean_mo_repro(subset(foliv, species == "SF" & sex == "F")); naoreprd = mean_mo_repro(subset(foliv, species == "NAO" & sex == "F"))
otreprd = mean_mo_repro(subset(insec, species == "OT" & sex == "F")); olreprd = mean_mo_repro(subset(insec, species == "OL" & sex == "F"))

# proportion of reproductive females by month and year
doreprdyr = mo_repro(subset(het, species == "DO" & sex == "F")); dmreprdyr = mo_repro(subset(het, species == "DM" & sex == "F")); pfreprdyr = mo_repro(subset(het, species == "PF" & sex == "F")); ppreprdyr = mo_repro(subset(het, species == "PP" & sex == "F")); pbreprdyr = mo_repro(subset(het, species == "PB" & sex == "F"));
pereprdyr = mo_repro(subset(cricet, species == "PE" & sex == "F")); pmreprdyr = mo_repro(subset(cricet, species == "PM" & sex == "F")); rmreprdyr = mo_repro(subset(cricet, species == "RM" & sex == "F"))
shreprdyr = mo_repro(subset(foliv, species == "SH" & sex == "F")); sfreprdyr = mo_repro(subset(foliv, species == "SF" & sex == "F")); naoreprdyr = mo_repro(subset(foliv, species == "NAO" & sex == "F"))
otreprdyr = mo_repro(subset(insec, species == "OT" & sex == "F")); olreprdyr = mo_repro(subset(insec, species == "OL" & sex == "F"))

#track the number of times females uniquely reproduce within years
doirep = indiv_repro(subset(het, species == "DO" & sex == "F")); dmirep = indiv_repro(subset(het, species == "DM" & sex == "F")); pfirep = indiv_repro(subset(het, species == "PF" & sex == "F")); ppirep = indiv_repro(subset(het, species == "PP" & sex == "F")); pbirep = indiv_repro(subset(het, species == "PB" & sex == "F"));
peirep = indiv_repro(subset(cricet, species == "PE" & sex == "F")); pmirep = indiv_repro(subset(cricet, species == "PM" & sex == "F")); rmirep = indiv_repro(subset(cricet, species == "RM" & sex == "F"))
shirep = indiv_repro(subset(foliv, species == "SH" & sex == "F")); sfirep = indiv_repro(subset(foliv, species == "SF" & sex == "F")); naoirep = indiv_repro(subset(foliv, species == "NAO" & sex == "F"))
otirep = indiv_repro(subset(insec, species == "OT" & sex == "F")); olirep = indiv_repro(subset(insec, species == "OL" & sex == "F"))

#proportion of years they were seen in
doyr = length(unique(het[het$species=="DO",]$yr))/10; dmyr = length(unique(het[het$species=="DM",]$yr))/10; pfyr = length(unique(het[het$species=="PF",]$yr))/10; ppyr = length(unique(het[het$species=="PP",]$yr))/10; pbyr = length(unique(het[het$species=="PB",]$yr))/10
peyr = length(unique(cricet[cricet$species=="PE",]$yr))/10; pmyr = length(unique(cricet[cricet$species=="PM",]$yr))/10; rmyr = length(unique(cricet[cricet$species=="RM",]$yr))/10
shyr = length(unique(foliv[foliv$species=="SH",]$yr))/10; sfyr = length(unique(foliv[foliv$species=="SF",]$yr))/10; naoyr = length(unique(foliv[foliv$species=="NAO",]$yr))/10
otyr = length(unique(insec[insec$species=="OT",]$yr))/10; olyr = length(unique(insec[insec$species=="OL",]$yr))/10

# average number of months they were seen in during years in which they were present
domo = mean_win_yr_occ(subset(het, species == "DO")); dmmo = mean_win_yr_occ(subset(het, species == "DM")); pfmo = mean_win_yr_occ(subset(het, species == "PF")); ppmo = mean_win_yr_occ(subset(het, species == "PP")); pbmo = mean_win_yr_occ(subset(het, species == "PB"))
pemo = mean_win_yr_occ(subset(cricet, species == "PE")); pmmo = mean_win_yr_occ(subset(cricet, species == "PM")); rmmo = mean_win_yr_occ(subset(cricet, species == "RM"))
shmo = mean_win_yr_occ(subset(foliv, species == "SH")); sfmo = mean_win_yr_occ(subset(foliv, species == "SF")); naomo = mean_win_yr_occ(subset(foliv, species == "NAO"))
otmo = mean_win_yr_occ(subset(insec, species == "OT")); olmo = mean_win_yr_occ(subset(insec, species == "OL"))

#mean abundance within all years 
doabun = allyrs_abun(subset(het, species == "DO")); dmabun = allyrs_abun(subset(het, species == "DM")); pfabun = allyrs_abun(subset(het, species == "PF")); ppabun = allyrs_abun(subset(het, species == "PP")); pbabun = allyrs_abun(subset(het, species == "PB"))
peabun = allyrs_abun(subset(cricet, species == "PE")); pmabun = allyrs_abun(subset(cricet, species == "PM")); rmabun = allyrs_abun(subset(cricet, species == "RM"))
shabun = allyrs_abun(subset(foliv, species == "SH")); sfabun = allyrs_abun(subset(foliv, species == "SF")); naoabun = allyrs_abun(subset(foliv, species == "NAO"))
otabun = allyrs_abun(subset(insec, species == "OT")); olabun = allyrs_abun(subset(insec, species == "OL"))

abuns = cbind(doabun, dmabun, pbabun, ppabun, pfabun, peabun, pmabun, rmabun, shabun, sfabun, naoabun, otabun, olabun)

#mean abundances within all years on control plots only
doabun = allyrs_abun(subset(het, species == "DO" & plot %in% c(1,2,4,8,9,11,12,14,17,22))); dmabun = allyrs_abun(subset(het, species == "DM"& plot %in% c(1,2,4,8,9,11,12,14,17,22))); pfabun = allyrs_abun(subset(het, species == "PF"& plot %in% c(1,2,4,8,9,11,12,14,17,22))); ppabun = allyrs_abun(subset(het, species == "PP"& plot %in% c(1,2,4,8,9,11,12,14,17,22))); pbabun = allyrs_abun(subset(het, species == "PB"& plot %in% c(1,2,4,8,9,11,12,14,17,22)))
peabun = allyrs_abun(subset(cricet, species == "PE"& plot %in% c(1,2,4,8,9,11,12,14,17,22))); pmabun = allyrs_abun(subset(cricet, species == "PM"& plot %in% c(1,2,4,8,9,11,12,14,17,22))); rmabun = allyrs_abun(subset(cricet, species == "RM"& plot %in% c(1,2,4,8,9,11,12,14,17,22)))
shabun = allyrs_abun(subset(foliv, species == "SH"& plot %in% c(1,2,4,8,9,11,12,14,17,22))); sfabun = allyrs_abun(subset(foliv, species == "SF"& plot %in% c(1,2,4,8,9,11,12,14,17,22))); naoabun = allyrs_abun(subset(foliv, species == "NAO"& plot %in% c(1,2,4,8,9,11,12,14,17,22)))
otabun = allyrs_abun(subset(insec, species == "OT"& plot %in% c(1,2,4,8,9,11,12,14,17,22))); olabun = allyrs_abun(subset(insec, species == "OL"& plot %in% c(1,2,4,8,9,11,12,14,17,22)))

conabuns = cbind(doabun, dmabun, pbabun, ppabun, pfabun, peabun, pmabun, rmabun, shabun, sfabun, naoabun, otabun, olabun)

#---------------------------------------------------------------------------------
#          calculate movement distances, multi-state capture histories
#---------------------------------------------------------------------------------

# get a vector unique tags, then get a vector of distances moved for all recaptured individuals, by SPECIES
  #heteromyids
    dmtags = unique(het[het$species == "DM",]$tag); dotags = unique(het[het$species == "DO",]$tag); pbtags = unique(het[het$species == "PB",]$tag); pptags = unique(het[het$species == "PP",]$tag); pftags = unique(het[het$species == "PF",]$tag)
    dmmeters = distance_moved(het[het$species == "DM",], dmtags); dometers = distance_moved(het[het$species == "DO",], dotags); pbmeters = distance_moved(het[het$species == "PB",], pbtags); ppmeters = distance_moved(het[het$species == "PP",], pptags); pfmeters = distance_moved(het[het$species == "PF",], pftags)
#cricetids
    petags = unique(cricet[cricet$species == "PE",]$tag); pmtags = unique(cricet[cricet$species == "PM",]$tag); rmtags = unique(cricet[cricet$species == "RM",]$tag)
    pemeters = distance_moved(cricet[cricet$species == "PE",], petags); pmmeters = distance_moved(cricet[cricet$species == "PM",], pmtags); rmmeters = distance_moved(cricet[cricet$species == "RM",], rmtags)
#folivores
    shtags = unique(foliv[foliv$species == "SH",]$tag); sftags = unique(foliv[foliv$species == "SF",]$tag); naotags = unique(foliv[foliv$species == "NAO",]$tag)
    shmeters = distance_moved(foliv[foliv$species == "SH",], shtags); sfmeters = distance_moved(foliv[foliv$species == "SF",], sftags); naometers = distance_moved(foliv[foliv$species == "NAO",], naotags)
#insectivores
    oltags = unique(insec[insec$species == "OL",]$tag); ottags = unique(insec[insec$species == "OT",]$tag)
    olmeters = distance_moved(insec[insec$species == "OL",], oltags); otmeters = distance_moved(insec[insec$species == "OT",], ottags)

# concatenate distance vectors for recaptured individuals by GUILD
Hgran = c(dmmeters, dometers, pbmeters, ppmeters, pfmeters)
Cgran = c(pemeters, pmmeters, rmmeters)
foli = c(shmeters, sfmeters) #separate NAO because they use different strategy - MIDDENS
insectiv = c(otmeters, olmeters)

corehet = c(dmmeters, dometers, pbmeters, ppmeters)

all = c(Hgran, Cgran, foli, naometers, insectiv)

# find breakpoints to use in MARK data structure for future analyses
# data reasonably well fits a lognormal distribution (eyeball and J. Powell)
# breakpoint = mean(logdata) + sd(logdata) of all the distances traveled by recaptured individuals    
  # using log1p, and back transforming using expm1 should solve the problem of lots of zeros 
corehet_brkpt = expm1(mean(log1p(corehet)) + sd(log1p(corehet)))
Hgran_brkpt = expm1(mean(log1p(Hgran)) + sd(log1p(Hgran)))
Cgran_brkpt = expm1(mean(log1p(Cgran)) + sd(log1p(Cgran)))
foli_brkpt = expm1(mean(log1p(foli)) + sd(log1p(foli)))
nao_brkpt = expm1(mean(log1p(naometers)) + sd(log1p(naometers)))
ins_brkpt = expm1(mean(log1p(insectiv)) + sd(log1p(insectiv)))

allbrk = expm1(mean(log1p(all)) + sd(log1p(all)))

# Get MARK capture histories
## add unique breakpoints for each species based on histogram data of movement
periods = c(261:380)
exclosures = c(5, 7, 10, 16, 23, 24)
krat_excl = c(5, 7, 10, 16, 23, 24, 3, 6, 13, 15, 18, 19, 20, 21)
DO_MARK = noplacelikehome(subset(het, species == "DO"), periods, krat_excl, corehet_brkpt)
DM_MARK = noplacelikehome(subset(het, species == "DM"), periods, krat_excl, corehet_brkpt)
PB_MARK = noplacelikehome(subset(het, species == "PB"), periods, exclosures, corehet_brkpt)
PP_MARK = noplacelikehome(subset(het, species == "PP"), periods, exclosures, corehet_brkpt)
PF_MARK = noplacelikehome(subset(het, species == "PF"), periods, exclosures, corehet_brkpt)
                          
PE_MARK = noplacelikehome(subset(cricet, species == "PE"), periods, exclosures, corehet_brkpt) 
PM_MARK = noplacelikehome(subset(cricet, species == "PM"), periods, exclosures, corehet_brkpt)
RM_MARK = noplacelikehome(subset(cricet, species == "RM"), periods, exclosures, corehet_brkpt)
                          
SH_MARK = noplacelikehome(subset(foliv, species == "SH"), periods, exclosures, corehet_brkpt)
SF_MARK = noplacelikehome(subset(foliv, species == "SF"), periods, exclosures, corehet_brkpt)
NAO_MARK = noplacelikehome(subset(foliv, species == "NAO"), periods, exclosures, corehet_brkpt)
                          
OT_MARK = noplacelikehome(subset(insec, species == "OT"), periods, exclosures, corehet_brkpt)
OL_MARK = noplacelikehome(subset(insec, species == "OL"), periods, exclosures, corehet_brkpt)

#---------------------------------------------------------------------------------
#          write files to folder for later analysis using RMark - .inp required
#---------------------------------------------------------------------------------

#write files to local folder
# write.table(DO_MARK, file = "mark_datafiles//do_mark.inp", row.names = F, col.names = F, quote = F)
# write.table(DM_MARK, file = "mark_datafiles//dm_mark.inp", row.names = F, col.names = F, quote = F)
# write.table(PB_MARK, file = "mark_datafiles//pb_mark.inp", row.names = F, col.names = F, quote = F)
# write.table(PP_MARK, file = "mark_datafiles//pp_mark.inp", row.names = F, col.names = F, quote = F)
# write.table(PF_MARK, file = "mark_datafiles//pf_mark.inp", row.names = F, col.names = F, quote = F)
# write.table(PE_MARK, file = "mark_datafiles//pe_mark.inp", row.names = F, col.names = F, quote = F)
# write.table(PM_MARK, file = "mark_datafiles//pm_mark.inp", row.names = F, col.names = F, quote = F)
# write.table(RM_MARK, file = "mark_datafiles//rm_mark.inp", row.names = F, col.names = F, quote = F)
# write.table(SH_MARK, file = "mark_datafiles//sh_mark.inp", row.names = F, col.names = F, quote = F)
# write.table(SF_MARK, file = "mark_datafiles//sf_mark.inp", row.names = F, col.names = F, quote = F)
# write.table(NAO_MARK, file = "mark_datafiles//nao_mark.inp", row.names = F, col.names = F, quote = F)
# write.table(OT_MARK, file = "mark_datafiles//ot_mark.inp", row.names = F, col.names = F, quote = F)
# write.table(OL_MARK, file = "mark_datafiles//ol_mark.inp", row.names = F, col.names = F, quote = F)

allspp = rbind(DO_MARK, DM_MARK, PB_MARK, PP_MARK, PF_MARK, PE_MARK, PM_MARK, RM_MARK, SH_MARK, SF_MARK, NAO_MARK, OT_MARK, OL_MARK)
write.table(allspp, file = "mark_datafiles//all_mark.inp", row.names = F, col.names = F, quote = F)


#---------------------------------------------------------------------------------
#          plot results
#---------------------------------------------------------------------------------


#------------------------------------------ FIGURE 2
pdf("Fig2_temporal_occ_spp_rank.pdf", 5, 5, pointsize = 10)
par(mfrow=c(1,1))

#average rank of each species 2000-2009 on control plots
#ordered: DM, DO, PB, PP, PF, PE, PM, RM, SH, SF, NA, OT, OL
conranks = c(2.9, 4, 2.6, 1, 7.875, 7.1, 9.8, 9.333, 9, 8.6, 8.5, 4.8, 10.5)

plot(conranks[5], pfyr, xlim = c(1,13), ylim = c(0,1), xlab = "species average ranked abundance", 
     ylab = "proportion of years present", pch = 19, col = "black", cex = 1.5, cex.axis = 1.5, cex.lab = 1.5)
  textxy(conranks[5], pfyr, "PF", cx = 1)
points(conranks[4],ppyr, pch = 19, col = "black", cex = 1.5)
  textxy(conranks[4],ppyr, "CP", cx = 1)
points(conranks[2], doyr, pch = 19, col = "black", cex = 1.5)
  textxy(conranks[2], doyr, "DO", cx = 1)
points(conranks[1], dmyr, pch = 19, col = "black", cex = 1.5)
  textxy(conranks[1], dmyr, "DM", cx = 1)
points(conranks[3], pbyr, pch = 19, col = "black", cex = 1.5)
  textxy(conranks[3], pbyr, "CB", cx = 1)
points(conranks[6], peyr, pch = 19, col = "gray30", cex = 1.5)
  textxy(conranks[6], peyr, "PE", cx = 1)
points(conranks[7], pmyr, pch = 19, col = "gray30", cex = 1.5)
  textxy(conranks[7], pmyr, "PM", cx = 1)
points(conranks[8], rmyr, pch = 19, col = "gray30", cex = 1.5)
  textxy(conranks[8], rmyr, "RM", cx = 1)
points(conranks[12], otyr, pch = 6, col = "gray30", cex = 1.5)
  textxy(conranks[12], otyr, "OT", cx = 1)
points(conranks[13], olyr, pch = 6, col = "gray30", cex = 1.5)
  textxy(conranks[13], olyr, "OL", cx = 1)
points(conranks[9], shyr, pch = 0, col = "gray30", cex = 1.5)
  textxy(conranks[9], shyr, "SH", cx = 1)
points(conranks[10], sfyr, pch = 0, col = "gray30", cex = 1.5)
  textxy(conranks[10], sfyr, "SF", cx = 1)
points(conranks[11], naoyr, pch = 0, col = "gray30", cex = 1.5)
  textxy(conranks[11], naoyr, "NA", cx = 1)

abline(v = 5, lty = 2, col = 'gray40', lwd = 1)
abline(h = 0.95, lty = 2, col = 'gray40', lwd = 1)

dev.off()


#------------------------------------------ FIGURE 3
pdf("Fig3_indiv_repro_trends.pdf", 6, 6, pointsize = 10)
par(mfrow=c(2,2))

#plot the relative abundance of females who represent each number of reproductive events per year
plot(NA, NA, xlim = c(0, 4), xaxp = c(0, 4, 4), xlab = "number reproductive events per year", 
     ylab = "proportion females", type = "b", pch = 19, ylim = c(0, 1), yaxp = c(0, 1, 4), bty = "n", 
     main = "heteromyid granivores", cex.axis = 1.5, cex.lab = 1.5)
points(c(0:3), table(doirep$num_reprod)/sum(table(doirep$num_reprod)), type = "b", pch = 19, col = "black")
points(c(0:4), table(dmirep$num_reprod)/sum(table(dmirep$num_reprod)), type = "b", pch = 15, col = "black")
points(c(0:4), table(pbirep$num_reprod)/sum(table(pbirep$num_reprod)), type = "b", pch = 17, col = "black")
points(c(0:4), table(ppirep$num_reprod)/sum(table(ppirep$num_reprod)), type = "b", pch = 25, col = "black")
points(c(0:2), table(pfirep$num_reprod)/sum(table(pfirep$num_reprod)), type = "b", pch = 24, col = "black")
legend('topright', c('DO', 'DM', 'CB', 'CP', 'PF'), bty = 'n', 
       col = 'black', pch = c(19,15,17,25,24), cex = 0.75)

plot(NA, NA, xlim = c(0, 4), xaxp = c(0, 4, 4), xlab = "number reproductive events per year", 
     ylab = "proportion females", type = "b", pch = 19, ylim = c(0, 1), yaxp = c(0, 1, 4), bty = "n", 
     main = "cricetid granivores", cex.axis = 1.5, cex.lab = 1.5)
points(c(0:2), table(peirep$num_reprod)/sum(table(peirep$num_reprod)), type = "b", pch = 19, col = "black")
points(c(0:2), table(pmirep$num_reprod)/sum(table(pmirep$num_reprod)), type = "b", pch = 15, col = "black")
points(c(0:1), table(rmirep$num_reprod)/sum(table(rmirep$num_reprod)), type = "b", pch = 17, col = "black")
legend('topright', c('PE', 'PM', 'RM'), bty = 'n', col = 'black', pch = c(19,15,17), cex = 0.75)

plot(NA, NA, xlim = c(0, 4), xaxp = c(0, 4, 4), xlab = "number reproductive events per year", 
     ylab = "proportion females", type = "b", pch = 19, ylim = c(0, 1), yaxp = c(0, 1, 4), bty = "n",
     main = "folivores", cex.axis = 1.5, cex.lab = 1.5)
points(c(0:1), table(shirep$num_reprod)/sum(table(shirep$num_reprod)), type = "b", pch = 19, col = "black")
points(c(0:1), table(sfirep$num_reprod)/sum(table(sfirep$num_reprod)), type = "b", pch = 15, col = "black")
points(c(0:3), table(naoirep$num_reprod)/sum(table(naoirep$num_reprod)), type = "b", pch = 17, col = "black")
legend('topright', c('SH', 'SF', 'NA'), bty = 'n', col = 'black', pch = c(19,15,17), cex = 0.75)

plot(NA, NA, xlim = c(0, 4), xaxp = c(0, 4, 4), xlab = "number reproductive events per year", bty = "n", 
     ylab = "proportion females", type = "b", pch = 19, ylim = c(0, 1), yaxp = c(0, 1, 4), 
     main = "carnivores", cex.axis = 1.5, cex.lab = 1.5)
points(c(0:3), table(otirep$num_reprod)/sum(table(otirep$num_reprod)), type = "b", pch = 19, col = "black")
points(c(0:1), table(olirep$num_reprod)/sum(table(olirep$num_reprod)), type = "b", pch = 15, col = "black")
legend('topright', c('OT', 'OL'), bty = 'n', col = 'black', pch = c(19,15), cex = 0.75)

dev.off()


#------------------------------------------ FIGURE 3
pdf("Fig4_species_movement_hist.pdf", 10, 8, pointsize = 10)
par(mfrow=c(5,3), mar=c(3,1.5,2,0.5), oma=c(1.5,2,1,1))

#plot histogram of all consecutive movement for rodents within a species 2000-2009
#create vector of breaks, incrementing by 6 meters (represents approx. 1 stake) since data are not actually continuous
v6 = seq(-3,500,6)
docount = hist(dometers, breaks = v6, col = 'gray60', xlim = c(0,500), ylim = c(0, 350), 
               xlab = "meters", main = 'DO - * - CORE')      
xline(corehet_brkpt, lwd = 2, col = "indianred")
xline(expm1(mean(log1p(dometers)) + sd(log1p(dometers))), lty = 3, lwd = 2) 

dmcount = hist(dmmeters, breaks = v6, col = 'gray60', xlim = c(0,500), ylim = c(0, 500), 
               xlab = "meters", main = 'DM - * - CORE')  
xline(corehet_brkpt, lwd = 2, col = "indianred")
xline(expm1(mean(log1p(dmmeters)) + sd(log1p(dmmeters))), lty = 3, lwd = 2) 

pbcount = hist(pbmeters, breaks = v6, col = 'gray60', xlim = c(0,500), ylim = c(0, 2000), 
               xlab = "meters", main = 'PB - * - CORE')   
xline(corehet_brkpt, lwd = 2, col = "indianred")
xline(expm1(mean(log1p(pbmeters)) + sd(log1p(pbmeters))), lty = 3, lwd = 2) 

ppcount = hist(ppmeters, breaks = v6, col = 'gray60', xlim = c(0,500), ylim = c(0, 600), 
               xlab = "meters", main = 'PP - * - CORE')      
xline(corehet_brkpt, lwd = 2, col = "indianred")
xline(expm1(mean(log1p(ppmeters)) + sd(log1p(ppmeters))), lty = 3, lwd = 2) 

otcount = hist(otmeters, breaks = v6, col = 'gray60', xlim = c(0,500), ylim = c(0,100), 
               xlab = "meters", main = 'OT - CORE')
xline(corehet_brkpt, lwd = 2, col = "indianred")
xline(expm1(mean(log1p(otmeters)) + sd(log1p(otmeters))), lty = 3, lwd = 2) 

pecount = hist(pemeters, breaks = v6, col = 'gray60', xlim = c(0,500), ylim = c(0,25), 
               xlab = "meters", main = 'PE - * - INT')
xline(corehet_brkpt, lwd = 2, col = "indianred")
xline(expm1(mean(log1p(pemeters)) + sd(log1p(pemeters))), lty = 3, lwd = 2) 

rmcount = hist(rmmeters, breaks = v6, col = 'gray60', xlim = c(0,500), ylim = c(0,15),
               xlab = "meters", main = 'RM - * - INT')
xline(corehet_brkpt, lwd = 2, col = "indianred")
xline(expm1(mean(log1p(rmmeters)) + sd(log1p(rmmeters))), lty = 3, lwd = 2) 

nacount = hist(naometers, breaks = v6, col = 'gray60', xlim = c(0,500), ylim = c(0,25), 
               xlab = "meters", main = 'NA - INT')
xline(corehet_brkpt, lwd = 2, col = "indianred")
xline(expm1(mean(log1p(naometers)) + sd(log1p(naometers))), lty = 3, lwd = 2) 

pfcount = hist(pfmeters, breaks = v6, col = 'gray60', xlim = c(0,500), ylim = c(0, 30), 
               xlab = "meters", main = 'PF - * - TRANS')  
xline(corehet_brkpt, lwd = 2, col = "indianred")
xline(expm1(mean(log1p(pfmeters)) + sd(log1p(pfmeters))), lty = 3, lwd = 2) 

pmcount = hist(pmmeters, breaks = v6, col = 'gray60', xlim = c(0,500), ylim = c(0,15),
               xlab = "meters", main = 'PM - * - TRANS')
xline(corehet_brkpt, lwd = 2, col = "indianred")
xline(expm1(mean(log1p(pmmeters)) + sd(log1p(pmmeters))), lty = 3, lwd = 2) 

shcount = hist(shmeters, breaks = v6, col = 'gray60', xlim = c(0,500), ylim = c(0,15), 
               xlab = "meters", main = 'SH - TRANS')
xline(corehet_brkpt, lwd = 2, col = "indianred")
xline(expm1(mean(log1p(shmeters)) + sd(log1p(shmeters))), lty = 3, lwd = 2) 

sfcount = hist(sfmeters, breaks = v6, col = 'gray60', xlim = c(0,500), ylim = c(0,15), 
               xlab = "meters", main = 'SF - TRANS')
xline(corehet_brkpt, lwd = 2, col = "indianred")
xline(expm1(mean(log1p(sfmeters)) + sd(log1p(sfmeters))), lty = 3, lwd = 2) 

olcount = hist(olmeters, breaks = v6, col = 'gray60', xlim = c(0,500), ylim = c(0,15), 
               xlab = "meters", main = 'OL - TRANS')
xline(corehet_brkpt, lwd = 2, col = "indianred")
xline(expm1(mean(log1p(olmeters)) + sd(log1p(olmeters))), lty = 3, lwd = 2) 

dev.off()


#--------------------------------------------------------------
#           FIGURES NOT INCLUDED IN THE MAIN PAPER
#--------------------------------------------------------------
pdf("Fig1_guild_movement_hist.pdf", 8, 5, pointsize = 10)
par(mfrow=c(3,2), mar=c(3,1.5,2,0.5), oma=c(1.5,2,1,1))

#plot histogram of all consecutive movement for rodents within a species 2000-2009
#create vector of breaks, incrementing by 6 meters (represents approx. 1 stake) since data are not actually continuous
v6 = seq(-3,500,6)
Hgcount = hist(Hgran, breaks = v6, col = 'gray60', xlim = c(0,500), ylim = c(0, 3500), 
               xlab = "meters", main = 'Heteromyids - PF, PP, PB, DO, DM')      
            xline(Hgran_brkpt, lwd = 2, col = "indianred")
Cgcount = hist(Cgran, breaks = v6, col = 'gray60', xlim = c(0,500), ylim = c(0,40), 
               xlab = "meters", main = 'Cricetids - PE, PM, RM')
            xline(Cgran_brkpt, lwd = 2, col = "indianred")
focount = hist(foli, breaks = v6, col = 'gray60', xlim = c(0,500), ylim = c(0,25),
               xlab = "meters", main = 'folivores - SH, SF')
            xline(foli_brkpt, lwd = 2, col = "indianred")
nacount = hist(naometers, breaks = v6, col = 'gray60', xlim = c(0,500), ylim = c(0,25), 
               xlab = "meters", main = 'neotoma - NA')
            xline(nao_brkpt, lwd = 2, col = "indianred")
incount = hist(insectiv, breaks = v6, col = 'gray60', xlim = c(0,500), ylim = c(0,120), 
               xlab = "meters", main = 'insectivores - OT, OL')
            xline(ins_brkpt, lwd = 2, col = "indianred")
dev.off()

#NOTES: Hgran has no obvious breakpoint (everything stays "home"), minima at 54 (distance between plots?)
#       Cgran has a minima at 48 m (distance between plots?)
#       Foli has a minima at 48 m (artefact of about the distance to hit between plots?)
#       Na has a minima at 24, closer to where I would guess Hgran to be (middens, never really leaves?)
#       Insec has a minima at 60 m, greater than the distance across a plot (active hunters)


#------------------------------------------------------------
pdf("Fig2_avg_temporal_occ.pdf", 5, 5, pointsize = 10)
par(mfrow=c(1,1))

# Make an occupancy plot for 2000-2009 (similar to Morgan) 
# plot temporal occupancy - for month and year 
plot(pfyr, pfmo, xlim = c(0,1), ylim = c(0,1), xlab = "across-year occupancy 2000-2009", ylab = "within-year occupancy", pch = 19, col = "hotpink")
    textxy(pfyr, pfmo, "PF")
  points(ppyr, ppmo, pch = 19, col = "hotpink")
    textxy(ppyr, ppmo, "PP")
  points(pbyr, pbmo, pch = 19, col = "hotpink")
    textxy(pbyr, pbmo, "PB")
  points(doyr, domo, pch = 19, col = "hotpink")
    textxy(doyr, domo, "DO")
  points(dmyr, dmmo, pch = 19, col = "hotpink")
    textxy(dmyr, dmmo, "DM")
  points(peyr, pemo, pch = 19, col = "indianred4")
    textxy(peyr, pemo, "PE")
  points(pmyr, pmmo, pch = 19. col = "indianred4")
    textxy(pmyr, pmmo, "PM")
  points(rmyr, rmmo, pch = 19, col = "indianred4")
    textxy(rmyr, rmmo, "RM")
  points(otyr, otmo, pch = 17, col = "cadetblue4")
    textxy(otyr, otmo, "OT")
  points(olyr, olmo, pch = 17, col = "cadetblue4")
    textxy(olyr, olmo, "OL")
  points(shyr, shmo, pch = 17, col = "chartreuse3")
    textxy(shyr, shmo, "SH")
  points(sfyr, sfmo, pch = 17, col = "chartreuse3")
    textxy(sfyr, sfmo, "SF")
  points(naoyr, naomo, pch = 17, col = "chartreuse3")
    textxy(naoyr, naomo, "NA")
  abline(v = 0.5, lty = 2, col = 'gray40', lwd = 1)
  abline(h = 0.5, lty = 2, col = 'gray40', lwd = 1)
dev.off()

#---------------------------------------------
pdf("Fig3_avg_prop_reprodfemales.pdf", 7, 10, paper = "letter", pointsize = 10)
par(mfrow=c(5,3))

#plot mean fecundity by month for each species
plot(c(1:12), doreprd, type = "l", xlim = c(1,12), ylim = c(0,1), pch = 19, xlab = "month", 
     ylab = "proprotion reproductive fem.", bty = "n", main = "DO - krat")
abline(h = 0.5, lty = 2, col = 'gray40', lwd = 1)

plot(c(1:12), dmreprd, type = "l", xlim = c(1,12), ylim = c(0,1), pch = 19, xlab = "month", 
     ylab = "proprotion reproductive fem.", bty = "n", main = "DM - krat")
abline(h = 0.5, lty = 2, col = 'gray40', lwd = 1)

plot(c(1:12), pbreprd, type = "l", xlim = c(1,12), ylim = c(0,1), pch = 19, xlab = "month", 
     ylab = "proprotion reproductive fem.", bty = "n", main = "PB - pocket mouse")
abline(h = 0.5, lty = 2, col = 'gray40', lwd = 1)

plot(c(1:12), ppreprd, type = "l", xlim = c(1,12), ylim = c(0,1), pch = 19, xlab = "month", 
     ylab = "proprotion reproductive fem.", bty = "n", main = "PP - pocket mouse")
abline(h = 0.5, lty = 2, col = 'gray40', lwd = 1)

plot(c(1:12), pfreprd, type = "l", xlim = c(1,12), ylim = c(0,1), pch = 19, xlab = "month", 
     ylab = "proprotion reproductive fem.", bty = "n", main = "PF - pocket mouse")
abline(h = 0.5, lty = 2, col = 'gray40', lwd = 1)

plot(c(1:12), pereprd, type = "l", xlim = c(1,12), ylim = c(0,1), pch = 19, xlab = "month", 
     ylab = "proprotion reproductive fem.", bty = "n", main = "PE - cactus mouse")
abline(h = 0.5, lty = 2, col = 'gray40', lwd = 1)

plot(c(1:12), pmreprd, type = "l", xlim = c(1,12), ylim = c(0,1), pch = 19, xlab = "month", 
     ylab = "proprotion reproductive fem.", bty = "n", main = "PM - deer mouse")
abline(h = 0.5, lty = 2, col = 'gray40', lwd = 1)

plot(c(1:12), rmreprd, type = "l", xlim = c(1,12), ylim = c(0,1), pch = 19, xlab = "month", 
     ylab = "proprotion reproductive fem.", bty = "n", main = "RM - harvest mouse")
abline(h = 0.5, lty = 2, col = 'gray40', lwd = 1)

plot(c(1:12), shreprd, type = "l", xlim = c(1,12), ylim = c(0,1), pch = 19, xlab = "month", 
     ylab = "proprotion reproductive fem.", bty = "n", main = "SH - cotton rat")
abline(h = 0.5, lty = 2, col = 'gray40', lwd = 1)

plot(c(1:12), sfreprd, type = "l", xlim = c(1,12), ylim = c(0,1), pch = 19, xlab = "month", 
     ylab = "proprotion reproductive fem.", bty = "n", main = "SF - cotton rat")
abline(h = 0.5, lty = 2, col = 'gray40', lwd = 1)

plot(c(1:12), naoreprd, type = "l", xlim = c(1,12), ylim = c(0,1), pch = 19, xlab = "month", 
     ylab = "proprotion reproductive fem.", bty = "n", main = "NA - pack rat")
abline(h = 0.5, lty = 2, col = 'gray40', lwd = 1)

plot(c(1:12), otreprd, type = "l", xlim = c(1,12), ylim = c(0,1), pch = 19, xlab = "month", 
     ylab = "proprotion reproductive fem.", bty = "n", main = "OT - grasshopper mouse")
abline(h = 0.5, lty = 2, col = 'gray40', lwd = 1)

plot(c(1:12), olreprd, type = "l", xlim = c(1,12), ylim = c(0,1), pch = 19, xlab = "month", 
     ylab = "proprotion reproductive fem.", bty = "n", main = "OL - grasshopper mouse")
abline(h = 0.5, lty = 2, col = 'gray40', lwd = 1)

dev.off()

#------------------------------------------
pdf("Fig4_abun_across_yrs.pdf", 8, 5, pointsize = 10)
par(mfrow=c(2,2), mar=c(2,2,2,0.5), oma=c(1.5,2,1,1))


years = c(2000:2009)

plot(years, doabun, type = "l", bty = "n", ylab = "", xlim = c(2000, 2010), 
     ylim = c(0, 600), xaxp = c(2000, 2009, 9))
  points(years, dmabun, type = "l")
  points(years, pbabun, type = "l")
  points(years, ppabun, type = "l")
  points(years, pfabun, type = "l")
    mtext("Heteromyidae granivores", side = 3, adj = 0.15, line = -1)
    mtext("abundance", side = 2, line = 2)

plot(years, peabun, type = "l", bty = "n", ylab = "", xlim = c(2000, 2010), 
     ylim = c(0, 600), xaxp = c(2000, 2009, 9))
points(years, pmabun, type = "l")
points(years, rmabun, type = "l")
    mtext("Cricetidae granivores", side = 3, adj = 0.15, line = -1)

plot(years, shabun, type = "l", bty = "n", ylab = "", xlim = c(2000, 2010), 
     ylim = c(0, 600), xaxp = c(2000, 2009, 9))
points(years, sfabun, type = "l")
points(years, naoabun, type = "l")
    mtext("Folivores", side = 3, adj = 0.15, line = -1)
    mtext("abundance", side = 2, line = 2)
    mtext("years", side = 1, line = 2)

plot(years, otabun, type = "l", bty = "n", ylab = "", xlim = c(2000, 2010), 
     ylim = c(0, 600), xaxp = c(2000, 2009, 9))
points(years, olabun, type = "l")
    mtext("Insectivores", side = 3, adj = 0.15, line = -1)
    mtext("years", side = 1, line = 2)

dev.off()


#------------------------------------------
pdf("Fig5_rank_abundance.pdf", 7, 10, paper = "letter", pointsize = 10)
par(mfrow=c(4,3))

years = c(2000:2009)
ranks = c(1:13)

#use data from control plots only
for (row in 1:nrow(conabuns)){
  yrdat = conabuns[row,]
  reldat = sort(yrdat/sum(yrdat), decreasing = TRUE)
    
  nonzero = reldat[reldat>0]
    labels = strtrim(as.character(names(nonzero)),2)
    print(labels)
  
  plot(ranks, reldat, type = "b", pch = 19, ylim = c(0,0.6), 
       xlab = "Rank", ylab = "Relative Abundance", bty = "n", xaxp = c(1, 13, 12))
  mtext(years[row], side = 3)
  textxy(c(1:length(nonzero)), nonzero, labs = labels, cx = 0.75)
}

dev.off()


#------------------------------------------
pdf("Fig7_mo-yr_repro_trends.pdf", 7, 10, paper = "letter", pointsize = 10)
par(mfrow=c(4,3))

alldat = rbind(doreprdyr, dmreprdyr, pbreprdyr, ppreprdyr, pfreprdyr, pereprdyr, pmreprdyr, rmreprdyr,
               shreprdyr, sfreprdyr, naoreprdyr, otreprdyr, olreprdyr)

years = c(2000:2009)
spp = c("DO", 'DM', 'PB', 'PP', 'PF', 'PE', 'PM', 'RM', 'SH', 'NAO', 'OT', 'OL') #SF (fit all on a page, hardly reproduce anyway)

#plot proportion fecundity by month and year for each species
for (y in 1:length (years)){
  
  for (s in 1:length (spp)) {
    data = subset(alldat, year == years[y] & species == spp[s])
    
    if(nrow(data) > 0){
  
    plot(c(1:12), data$proprepro, type = "l", xlim = c(1,12), ylim = c(0,1), pch = 19, xlab = "month",  
       ylab = "prop. reproductive F", bty = "n", main =  paste(spp[s], years[y], sep = " "))
    abline(h = 0.5, lty = 2, col = 'gray40', lwd = 1)
    
}
  else {
    plot(NA, NA, type = "l", xlim = c(1,12), ylim = c(0,1), pch = 19, xlab = "month",  
         ylab = "prop. reproductive F", bty = "n",  main = paste( "No", spp[s], years[y], sep = " "))
  }}}

dev.off()


#------------------------------------------
pdf("Fig8_histogram_masses.pdf", 6, 6, pointsize = 10)
par(mfrow=c(2,2))

v10 = seq(0,300,10)
hist(het$wgt, col = 'gray60', xlim = c(0,80), ylim = c(0, 3000), xlab = "grams", 
     main = "Heteromyid mass")

hist(cricet$wgt, col = 'gray60', xlim = c(0,40), ylim = c(0, 500), xlab = "grams", 
     main = "Cricetid granivore mass")

hist(foliv$wgt, col = 'gray60', xlim = c(0,300), ylim = c(0, 125), xlab = "grams", 
     main = "folivore mass")

hist(insec$wgt, col = 'gray60', xlim = c(0,50), ylim = c(0, 800), xlab = "grams", 
     main = "carnivore mass")

dev.off()

#------------------------------------------
#proportion of years they were seen in
doyr = length(unique(allrats[allrats$species=="DO",]$yr))/35; dmyr = length(unique(allrats[allrats$species=="DM",]$yr))/35; pfyr = length(unique(allrats[allrats$species=="PF",]$yr))/35; ppyr = length(unique(allrats[allrats$species=="PP",]$yr))/35; pbyr = length(unique(allrats[allrats$species=="PB",]$yr))/35
peyr = length(unique(allrats[allrats$species=="PE",]$yr))/35; pmyr = length(unique(allrats[allrats$species=="PM",]$yr))/35; rmyr = length(unique(allrats[allrats$species=="RM",]$yr))/35
shyr = length(unique(allrats[allrats$species=="SH",]$yr))/35; sfyr = length(unique(allrats[allrats$species=="SF",]$yr))/35; naoyr = length(unique(allrats[allrats$species=="NAO",]$yr))/35
otyr = length(unique(allrats[allrats$species=="OT",]$yr))/35; olyr = length(unique(allrats[allrats$species=="OL",]$yr))/35

# average number of months they were seen in during years in which they were present
domo = mean_win_yr_alldat(subset(allrats, species == "DO")); dmmo = mean_win_yr_alldat(subset(allrats, species == "DM")); pfmo = mean_win_yr_alldat(subset(allrats, species == "PF")); ppmo = mean_win_yr_alldat(subset(allrats, species == "PP")); pbmo = mean_win_yr_alldat(subset(allrats, species == "PB"))
pemo = mean_win_yr_alldat(subset(allrats, species == "PE")); pmmo = mean_win_yr_alldat(subset(allrats, species == "PM")); rmmo = mean_win_yr_alldat(subset(allrats, species == "RM"))
shmo = mean_win_yr_alldat(subset(allrats, species == "SH")); sfmo = mean_win_yr_alldat(subset(allrats, species == "SF")); naomo = mean_win_yr_alldat(subset(allrats, species == "NAO"))
otmo = mean_win_yr_alldat(subset(allrats, species == "OT")); olmo = mean_win_yr_alldat(subset(allrats, species == "OL"))

#mean abundance within all years 
doabun = allyrs_abun(subset(allrats, species == "DO")); dmabun = allyrs_abun(subset(het, species == "DM")); pfabun = allyrs_abun(subset(het, species == "PF")); ppabun = allyrs_abun(subset(het, species == "PP")); pbabun = allyrs_abun(subset(het, species == "PB"))
peabun = allyrs_abun(subset(cricet, species == "PE")); pmabun = allyrs_abun(subset(cricet, species == "PM")); rmabun = allyrs_abun(subset(cricet, species == "RM"))
shabun = allyrs_abun(subset(foliv, species == "SH")); sfabun = allyrs_abun(subset(foliv, species == "SF")); naoabun = allyrs_abun(subset(foliv, species == "NAO"))
otabun = allyrs_abun(subset(insec, species == "OT")); olabun = allyrs_abun(subset(insec, species == "OL"))

# Make an occupancy plot for 1978-2012 (similar to Morgan) 
# plot temporal occupancy - for month and year 
pdf("Fig9_avg_temporal_occ_allyears.pdf", 5, 5, pointsize = 10)
par(mfrow=c(1,1))

plot(pfyr, pfmo, xlim = c(0,1), ylim = c(0,1), xlab = "across-year occupancy 1978-2012", ylab = "within-year occupancy", pch = 19, col = "hotpink")
textxy(pfyr, pfmo, "PF")
points(ppyr, ppmo, pch = 19, col = "hotpink")
textxy(ppyr, ppmo, "PP")
points(pbyr, pbmo, pch = 19, col = "hotpink")
textxy(pbyr, pbmo, "PB")
points(doyr, domo, pch = 19, col = "hotpink")
textxy(doyr, domo, "DO")
points(dmyr, dmmo, pch = 19, col = "hotpink")
textxy(dmyr, dmmo, "DM")
points(peyr, pemo, pch = 19, col = "indianred4")
textxy(peyr, pemo, "PE")
points(pmyr, pmmo, pch = 19, col = "indianred4")
textxy(pmyr, pmmo, "PM")
points(rmyr, rmmo, pch = 19, col = "indianred4")
textxy(rmyr, rmmo, "RM")
points(otyr, otmo, pch = 19, col = "cadetblue4")
textxy(otyr, otmo, "OT")
points(olyr, olmo, pch = 17, col = "cadetblue4")
textxy(olyr, olmo, "OL")
points(shyr, shmo, pch = 17, col = "chartreuse3")
textxy(shyr, shmo, "SH")
points(sfyr, sfmo, pch = 17, col = "chartreuse3")
textxy(sfyr, sfmo, "SH")
points(naoyr, naomo, pch = 17, col = "chartreuse3")
textxy(naoyr, naomo, "NA")

abline(v = 0.5, lty = 2, col = 'gray40', lwd = 1)
abline(h = 0.5, lty = 2, col = 'gray40', lwd = 1)

dev.off()


#------------------------------------------------------
pdf("Fig12_avg_repro_trends_wpts.pdf", 7, 10, paper = "letter", pointsize = 10)
par(mfrow=c(5,3))

alldat = rbind(doreprdyr, dmreprdyr, pbreprdyr, ppreprdyr, pfreprdyr, pereprdyr, pmreprdyr, rmreprdyr,
               shreprdyr, sfreprdyr, naoreprdyr, otreprdyr, olreprdyr)

mos = c(1:12)
spp = c("DO", 'DM', 'PB', 'PP', 'PF', 'PE', 'PM', 'RM', 'SH', 'SF', 'NAO', 'OT', 'OL') 

#plot proportion fecundity by month and year for each species
for (s in 1:length (spp)) {
  avg_repro = c()
  z = c() # for determining color of points
  
  for (m in 1:length (mos)){
  
    data = subset(alldat, month == mos[m] & species == spp[s])
    mo_repro = as.numeric(data$proprepro)
    mo_repro = mo_repro[!is.na(mo_repro)] #remove NAs from data
    
    if (length(mo_repro > 0)){
      avg = mean(mo_repro)
      
      if (avg >= 0.5) {
        prop_mos = length(mo_repro[mo_repro > 0.5])/length(mo_repro)  #proportion of months in which the spp was seen that had repro > 50%     
        z = append(z, prop_mos)
      }
      else {
        z = append(z, NA)
      }
    }
      else {
      avg = NA }
      
    avg_repro = append(avg_repro, avg)   
  }
  
  plot(mos, avg_repro, type = "l", xlim = c(1,12), ylim = c(0,1), pch = 19, xlab = "month", 
         ylab = "proprotion reproductive fem.", bty = "n", main = spp[s])
    abline(h = 0.5, lty = 2, col = 'gray40', lwd = 1)
  
  points(mos, avg_repro, col = ifelse(z >=0.5, "indianred", "cadetblue"), pch = 19)
}

dev.off()


#-------------------------------------------------------------------
#          Print statments - descriptive info for the txt file
#-------------------------------------------------------------------

#close sink file
sink()
# 


# 
# ##------ histograms of the number of recaptures (some species are recaptured far less often), save in pdf file #FIXME
# pdf(file = "recaps_by_species.pdf", 11, 7.5)
# par(mfrow = c(2,4))
# 
# plot_recap_hist(PB, "PB"); plot_recap_hist(PP, "PP"); plot_recap_hist(PF, "PF"); plot_recap_hist(DM, "DM"); plot_recap_hist(DO, "DO")
# plot_recap_hist(PE, "PE"); plot_recap_hist(PM, "PM"); plot_recap_hist(RM, "RM")
# plot_recap_hist(NAO, "NAO"); plot_recap_hist(SF, "SF"); plot_recap_hist(SH, "SH")
# plot_recap_hist(OT, "OT"); plot_recap_hist(OL, "OL")
# 
# dev.off()

