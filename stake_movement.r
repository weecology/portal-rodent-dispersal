# Code for working with individual-level rodent data
# movement with stakes

library(calibrate)
library(fields)

#---------------------------------------------------------------------------------
#          setup - select wd, import data, source code,  file to collect results
#---------------------------------------------------------------------------------
#set working directory
wd = "C://Users//sarah//Documents//GitHub//portal-rodent-dispersal"
setwd(wd)
source("movement_fxns.R")

#import the data by guild
het = read.csv("rawdata//heteromyids_2000-2009.csv")   # DO, DM, PB, PP, PF
cricet = read.csv("rawdata//cricetids_2000-2009.csv")  # PE, PM, RM
foliv = read.csv("rawdata//folivores_2000-2009.csv")   # SH, SF, NA (as NAO)
insec = read.csv("rawdata//onychomys_2000-2009.csv")   # OT, OL

# output directed to rodent_results.txt in wd. output is appended
# to existing file. output also sent to terminal. 
sink("rodent_results.txt", append=TRUE, split=TRUE)

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

# find breakpoints to use in MARK data structure for future analyses
# data reasonably well fits a lognormal distribution (eyeball and J. Powell)
# breakpoint = mean(logdata) + sd(logdata) of all the distances traveled by recaptured individuals    
Hgran_brkpt = expm1(mean(log1p(Hgran)) + sd(log1p(Hgran)))
Cgran_brkpt = expm1(mean(log1p(Cgran)) + sd(log1p(Cgran)))
foli_brkpt = expm1(mean(log1p(foli)) + sd(log1p(foli)))
nao_brkpt = expm1(mean(log1p(naometers)) + sd(log1p(naometers)))
ins_brkpt = expm1(mean(log1p(insectiv)) + sd(log1p(insectiv)))

# Get MARK capture histories
## add unique breakpoints for each species based on histogram data of movement
periods = c(261:380)
exclosures = c(5, 7, 10, 16, 23, 24)
krat_excl = c(5, 7, 10, 16, 23, 24, 3, 6, 13, 15, 18, 19, 20, 21)
DO_MARK = noplacelikehome(subset(het, species == "DO"), periods, krat_excl, Hgran_brkpt)
DM_MARK = noplacelikehome(subset(het, species == "DM"), periods, krat_excl, Hgran_brkpt)
PB_MARK = noplacelikehome(subset(het, species == "PB"), periods, exclosures, Hgran_brkpt)
PP_MARK = noplacelikehome(subset(het, species == "PP"), periods, exclosures, Hgran_brkpt)
PF_MARK = noplacelikehome(subset(het, species == "PF"), periods, exclosures, Hgran_brkpt)
                          
PE_MARK = noplacelikehome(subset(cricet, species == "PE"), periods, exclosures, Cgran_brkpt) 
PM_MARK = noplacelikehome(subset(cricet, species == "PM"), periods, exclosures, Cgran_brkpt)
RM_MARK = noplacelikehome(subset(cricet, species == "RM"), periods, exclosures, Cgran_brkpt)
                          
SH_MARK = noplacelikehome(subset(foliv, species == "SH"), periods, exclosures, foli_brkpt)
SF_MARK = noplacelikehome(subset(foliv, species == "SF"), periods, exclosures, foli_brkpt)

NAO_MARK = noplacelikehome(subset(foliv, species == "NAO"), periods, exclosures, nao_brkpt)
                          
OT_MARK = noplacelikehome(subset(insec, species == "OT"), periods, exclosures, ins_brkpt)
OL_MARK = noplacelikehome(subset(insec, species == "OL"), periods, exclosures, ins_brkpt)

#---------------------------------------------------------------------------------
#          write files to folder for later analysis using RMark - .inp required
#---------------------------------------------------------------------------------

#write files to local folder
write.table(DO_MARK, file = "mark_datafiles//do_mark.inp", row.names = F, col.names = F)
write.table(DM_MARK, file = "mark_datafiles//dm_mark.inp", row.names = F, col.names = F)
write.table(PB_MARK, file = "mark_datafiles//pb_mark.inp", row.names = F, col.names = F)
write.table(PP_MARK, file = "mark_datafiles//pp_mark.inp", row.names = F, col.names = F)
write.table(PF_MARK, file = "mark_datafiles//pf_mark.inp", row.names = F, col.names = F)
write.table(PE_MARK, file = "mark_datafiles//pe_mark.inp", row.names = F, col.names = F)
write.table(PM_MARK, file = "mark_datafiles//pm_mark.inp", row.names = F, col.names = F)
write.table(RM_MARK, file = "mark_datafiles//rm_mark.inp", row.names = F, col.names = F)
write.table(SH_MARK, file = "mark_datafiles//sh_mark.inp", row.names = F, col.names = F)
write.table(SF_MARK, file = "mark_datafiles//sf_mark.inp", row.names = F, col.names = F)
write.table(NAO_MARK, file = "mark_datafiles//nao_mark.inp", row.names = F, col.names = F)
write.table(OT_MARK, file = "mark_datafiles//ot_mark.inp", row.names = F, col.names = F)
write.table(OL_MARK, file = "mark_datafiles//ol_mark.inp", row.names = F, col.names = F)

#---------------------------------------------------------------------------------
#          plot results
#---------------------------------------------------------------------------------
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


#### Make an occupancy plot for 2000-2009 (similar to Morgan) 
#proportion of years they were seen in
doyr = length(unique(het[het$species=="DO",]$yr))/10; dmyr = length(unique(het[het$species=="DM",]$yr))/10; pfyr = length(unique(het[het$species=="PF",]$yr))/10; ppyr = length(unique(het[het$species=="PP",]$yr))/10; pbyr = length(unique(het[het$species=="PB",]$yr))/10
peyr = length(unique(cricet[cricet$species=="PE",]$yr))/10; pmyr = length(unique(cricet[cricet$species=="PM",]$yr))/10; rmyr = length(unique(cricet[cricet$species=="RM",]$yr))/10
shyr = length(unique(foliv[foliv$species=="SH",]$yr))/10; sfyr = length(unique(foliv[foliv$species=="SF",]$yr))/10; naoyr = length(unique(foliv[foliv$species=="NAO",]$yr))/10
otyr = length(unique(insec[insec$species=="OT",]$yr))/10; olyr = length(unique(insec[insec$species=="OL",]$yr))/10

# average number of months they were seen in during year sin which they were present
domo = mean_win_yr_occ(subset(het, species == "DO")); dmmo = mean_win_yr_occ(subset(het, species == "DM")); pfmo = mean_win_yr_occ(subset(het, species == "PF")); ppmo = mean_win_yr_occ(subset(het, species == "PP")); pbmo = mean_win_yr_occ(subset(het, species == "PB"))
pemo = mean_win_yr_occ(subset(cricet, species == "PE")); pmmo = mean_win_yr_occ(subset(cricet, species == "PM")); rmmo = mean_win_yr_occ(subset(cricet, species == "RM"))
shmo = mean_win_yr_occ(subset(foliv, species == "SH")); sfmo = mean_win_yr_occ(subset(foliv, species == "SF")); naomo = mean_win_yr_occ(subset(foliv, species == "NAO"))
otmo = mean_win_yr_occ(subset(insec, species == "OT")); olmo = mean_win_yr_occ(subset(insec, species == "OL"))
 
# plot temporal occupancy - for month and year 
pdf("Fig2_avg_temporal_occ.pdf", 5, 5, pointsize = 10)
par(mfrow=c(1,1))

plot(pfyr, pfmo, xlim = c(0,1), ylim = c(0,1), xlab = "across-year occupancy", ylab = "within-year occupancy", pch = 19, col = "hotpink")
    textxy(pfyr, pfmo, "PF")
  points(ppyr, ppmo, pch = 19, col = "hotpink")
    textxy(ppyr, ppmo, "PP")
  points(pbyr, pbmo, pch = 19, col = "hotpink")
    textxy(pbyr, pbmo, "PB")
  points(doyr, domo, pch = 19, col = "hotpink")
    textxy(doyr, domo, "DO")
  points(dmyr, dmmo, pch = 19, col = "hotpink")
    textxy(dmyr, dmmo, "DM")
  points(peyr, pemo, pch = 19)
    textxy(peyr, pemo, "Cgran")
  points(pmyr, pmmo, pch = 19)
    textxy(pmyr, pmmo, "Cgran")
  points(otyr, otmo, pch = '*', cex = 1.5)
    textxy(otyr, otmo, "insectiv")
  points(olyr, olmo, pch = '*', cex = 1.5)
    textxy(olyr, olmo, "insectiv")
  points(shyr, shmo, pch = 19)
    textxy(shyr, shmo, "foliv")
  points(sfyr, sfmo, pch = 19)
    textxy(sfyr, sfmo, "foliv")
  points(naoyr, naomo, pch = 19, col = "purple")
    textxy(naoyr, naomo, "NA")
  points(rmyr, rmmo, pch = 19)
    textxy(rmyr, rmmo, "Cgran")
  abline(v = 0.5, lty = 2, col = 'gray40', lwd = 1)
  abline(h = 0.5, lty = 2, col = 'gray40', lwd = 1)

dev.off()

# #mean abundance within all years #FIXME - CAN'T FIND FXN
# doabun = mean_allyrs_abun(subset(het, species == "DO")); dmabun = mean_allyrs_abun(DM); pfabun = mean_allyrs_abun(PF); ppabun = mean_allyrs_abun(PP); pbabun = mean_allyrs_abun(PB)
# peabun = mean_allyrs_abun(PE); pmabun = mean_allyrs_abun(PM); rmabun = mean_allyrs_abun(RM)
# shabun = mean_allyrs_abun(SH); sfabun = mean_allyrs_abun(SF); naoabun = mean_allyrs_abun(NAO)
# otabun = mean_allyrs_abun(OT); olabun = mean_allyrs_abun(OL)
# 

# ##------- plot number of individuals captured in each period for each group, save in pdf file #FIXME
# pdf(file="indivs_per_year.pdf",11,7.5)
# par(mfrow=c(3,2))
# 
# plot_freq_by_prd(subset(het, species == PB), "PB"); plot_freq_by_prd(PP, "PP"); plot_freq_by_prd(PF, "PF"); plot_freq_by_prd(DM, "DM"); plot_freq_by_prd(DO, "DO")
# plot_freq_by_prd(PE, "PE"); plot_freq_by_prd(PM, "PM"); plot_freq_by_prd(RM, "RM")
# plot_freq_by_prd(NAO, "NAO"); plot_freq_by_prd(SF, "SF"); plot_freq_by_prd(SH, "SH")
# plot_freq_by_prd(OT, "OT"); plot_freq_by_prd(OL, "OL")
# 
# dev.off()
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

