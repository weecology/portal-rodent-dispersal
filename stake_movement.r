# Code for working with individual-level rodent data
# movement with stakes

library(calibrate)

wd = "C://Users//sarah//Desktop//Dropbox//Active Research Projects//Rodent Movement"
setwd(wd)
source("C://Users//sarah//Documents//GitHub//portal-rodent-dispersal//movement_fxns.R")

het = read.csv("data//heteromyids_2000-2009.csv")   # DO, DM, PB, PP, PF
cricet = read.csv("data//cricetids_2000-2009.csv")  # PE, PM, RM
foliv = read.csv("data//folivores_2000-2009.csv")   # SH, SF, NA (as NAO)
insec = read.csv("data//onychomys_2000-2009.csv")   # OT, OL

# # # list of dataframes                                                 WORKS NOW, TO USE? OR NOT TO USE? 
# datLS = list(het, cricet, foliv, insec)
# names(datLS) = c('het','cricet', 'foliv', 'insec')
# for (df in seq_along(datLS)) datLS[[df]]$tag = as.character(datLS[[df]]$tag)


# change some cols from factor to character class
het$tag = as.character(het$tag); cricet$tag = as.character(cricet$tag); foliv$tag = as.character(foliv$tag); insec$tag = as.character(insec$tag)
het$species = as.character(het$species); cricet$species = as.character(cricet$species); foliv$species = as.character(foliv$species); insec$species = as.character(insec$species)
het$sex = as.character(het$sex); cricet$sex = as.character(cricet$sex); foliv$sex = as.character(foliv$sex); insec$sex = as.character(insec$sex)

# give untagged individuals a unique 7-number code
het = id_unknowns(het, 16); cricet = id_unknowns(cricet, 16); foliv = id_unknowns(foliv, 16); insec = id_unknowns(insec, 16)

# get rid of 'bad data'; deletes data that is not a pit tag, where sex is inconsistent or where species is inconsistent. 
het = subsetDat(het); cricet = subsetDat(cricet); foliv = subsetDat(foliv); insec = subsetDat(insec)


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
# breakpoint = mean + sd of all the distances traveled by recaptured individuals    #IS THERE A BETTER WAY!?
Hgran_brkpt = mean(Hgran) + sd(Hgran)
Cgran_brkpt = mean(Cgran) + sd(Cgran)
foli_brkpt = mean(foli) + sd(foli)
nao_brkpt = mean(naometers) + sd(naometers)
ins_brkpt = mean(insectiv) + sd(insectiv)

#plot histogram of all consecutive movement for rodents within a species 2000-2009
#create vector of breaks, incrementing by 6 meters (represents approx. 1 stake) since data are not actually continuous
v6 = seq(-3,500,6)
Hgcount = hist(Hgran, breaks = v6, col = 'gray60', xlim = c(0,500), ylim = c(0, 2500), main = 'Heteromyids - PF, PP, PB, DO, DM')      
  xline(Hgran_brkpt, lwd = 2, col = "indianred")
Cgcount = hist(Cgran, breaks = v6, col = 'gray60', xlim = c(0,500), ylim = c(0,20), main = 'Cricetids - PE, PM, RM')
  xline(Cgran_brkpt, lwd = 2, col = "indianred")
focount = hist(foli, breaks = v6, col = 'gray60', xlim = c(0,500), main = 'folivores - SH, SF')
  xline(foli_brkpt, lwd = 2, col = "indianred")
nacount = hist(naometers, breaks = v6, col = 'gray60', xlim = c(0,500), main = 'neotoma - NA')
  xline(nao_brkpt, lwd = 2, col = "indianred")
incount = hist(insectiv, breaks = v6, col = 'gray60', xlim = c(0,500), ylim = c(0,80), main = 'insectivores - OT, OL')
  xline(ins_brkpt, lwd = 2, col = "indianred")

# Get MARK capture histories
## add unique breakpoints for each species based on histogram data of movement
periods = c(261:380)
DO_MARK = noplacelikehome(subset(het, species == "DO"), periods, Hgran_brkpt)
DM_MARK = noplacelikehome(subset(het, species == "DM"), periods, Hgran_brkpt)
PB_MARK = noplacelikehome(subset(het, species == "PB"), periods, Hgran_brkpt)
PP_MARK = noplacelikehome(subset(het, species == "PP"), periods, Hgran_brkpt)
PF_MARK = noplacelikehome(subset(het, species == "PF"), periods, Hgran_brkpt)
                          
PE_MARK = noplacelikehome(subset(cricet, species == "PE"), periods, Cgran_brkpt) 
PM_MARK = noplacelikehome(subset(cricet, species == "PM"), periods, Cgran_brkpt)
RM_MARK = noplacelikehome(subset(cricet, species == "RM"), periods, Cgran_brkpt)
                          
SH_MARK = noplacelikehome(subset(foliv, species == "SH"), periods, foli_brkpt)
SF_MARK = noplacelikehome(subset(foliv, species == "SF"), periods, foli_brkpt)

NAO_MARK = noplacelikehome(subset(foliv, species == "NAO"), periods, nao_brkpt)
                          
OT_MARK = noplacelikehome(subset(insec, species == "OT"), periods, ins_brkpt)
OL_MARK = noplacelikehome(subset(insec, species == "OL"), periods, ins_brkpt)


############### PLOTTING

# plot density of movment by guild for 2000-2009
plot(density(Hgran), main = 'Portal movement by guild', xlab = 'meters', lwd = 2, col = 'hotpink2')
  lines(density(Cgran), col = 'deepskyblue3', lwd = 3, lty = 6)
  lines(density(naometers), col = 'indianred4', lwd = 4, lty = 3)
  lines(density(foli), col = 'mediumpurple4', lwd = 4, lty = 3)
  lines(density(insectiv), col = 'darkgreen', lwd = 2)
  legend('topright', c('Hgran', 'Neotoma', 'foliv', 'Cgran', 'insec'), bty = 'n', lty = c(1,3,6,3,1), lwd = 5, seg.len = 2,
         col = c('hotpink2', 'indianred4', 'mediumpurple4', 'deepskyblue3', 'darkgreen'))

#### Make an occupancy plot for 2000-2009 (similar to Morgan) #FIX ME - REFERS TO OLD DATAFRAMES
#proportion of years they were seen in
doyr = length(unique(DO$yr))/10; dmyr = length(unique(DM$yr))/10; pfyr = length(unique(PF$year))/10; ppyr = length(unique(PP$year))/10; pbyr = length(unique(PB$year))/10
peyr = length(unique(PE$yr))/10; pmyr = length(unique(PM$yr))/10; rmyr = length(unique(RM$yr))/10
shyr = length(unique(SH$yr))/10; sfyr = length(unique(SF$yr))/10; naoyr = length(unique(NAO$yr))/10
otyr = length(unique(OT$yr))/10; olyr = length(unique(OL$yr))/10

#proportion of within-year trapping periods they were seen in 
domo = mean_win_yr_occ(DO); dmmo = mean_win_yr_occ(DM); pfmo = mean_win_yr_occ(PF); ppmo = mean_win_yr_occ(PP); pbmo = mean_win_yr_occ(PB)
pemo = mean_win_yr_occ(PE); pmmo = mean_win_yr_occ(PM); rmmo = mean_win_yr_occ(RM)
shmo = mean_win_yr_occ(SH); sfmo = mean_win_yr_occ(SF); naomo = mean_win_yr_occ(NAO)
otmo = mean_win_yr_occ(OT); olmo = mean_win_yr_occ(OL)

#mean abundance within all years
doabun = mean_allyrs_abun(DO); dmabun = mean_allyrs_abun(DM); pfabun = mean_allyrs_abun(PF); ppabun = mean_allyrs_abun(PP); pbabun = mean_allyrs_abun(PB)
peabun = mean_allyrs_abun(PE); pmabun = mean_allyrs_abun(PM); rmabun = mean_allyrs_abun(RM)
shabun = mean_allyrs_abun(SH); sfabun = mean_allyrs_abun(SF); naoabun = mean_allyrs_abun(NAO)
otabun = mean_allyrs_abun(OT); olabun = mean_allyrs_abun(OL)


# plot temporal occupancy - for month and year
plot(pfyr, pfmo, xlim = c(0,1), ylim = c(0,1), xlab = "across-year occupancy", ylab = "within-year occupancy", pch = 19, col = "hotpink")
    #textxy(pfyr, pfmo, "PF")
  points(ppyr, ppmo, pch = 19, col = "hotpink")
    #textxy(ppyr, ppmo, "PP")
  points(pbyr, pbmo, pch = 19, col = "hotpink")
    #textxy(pbyr, pbmo, "PB")
  points(doyr, domo, pch = 19, col = "hotpink")
    #textxy(doyr, domo, "DO")
  points(dmyr, dmmo, pch = 19, col = "hotpink")
    #textxy(dmyr, dmmo, "DM")
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
  abline(v = 0.5, lty = 2, col = 'gray40', lwd = 2)
  abline(h = 0.5, lty = 2, col = 'gray40', lwd = 2)

##------- plot number of individuals captured in each period for each group, save in pdf file
pdf(file="indivs_per_year.pdf",11,7.5)
par(mfrow=c(3,2))

plot_freq_by_prd(PB, "PB"); plot_freq_by_prd(PP, "PP"); plot_freq_by_prd(PF, "PF"); plot_freq_by_prd(DM, "DM"); plot_freq_by_prd(DO, "DO")
plot_freq_by_prd(PE, "PE"); plot_freq_by_prd(PM, "PM"); plot_freq_by_prd(RM, "RM")
plot_freq_by_prd(NAO, "NAO"); plot_freq_by_prd(SF, "SF"); plot_freq_by_prd(SH, "SH")
plot_freq_by_prd(OT, "OT"); plot_freq_by_prd(OL, "OL")

dev.off()

##------ histograms of the number of recaptures (some species are recaptured far less often), save in pdf file
pdf(file = "recaps_by_species.pdf", 11, 7.5)
par(mfrow = c(2,4))

plot_recap_hist(PB, "PB"); plot_recap_hist(PP, "PP"); plot_recap_hist(PF, "PF"); plot_recap_hist(DM, "DM"); plot_recap_hist(DO, "DO")
plot_recap_hist(PE, "PE"); plot_recap_hist(PM, "PM"); plot_recap_hist(RM, "RM")
plot_recap_hist(NAO, "NAO"); plot_recap_hist(SF, "SF"); plot_recap_hist(SH, "SH")
plot_recap_hist(OT, "OT"); plot_recap_hist(OL, "OL")

dev.off()

######################### EXTRA STUFF, PROBABLY DON'T NEED.... ##############################


############## MOVING
# get list of indivs that moved plots or treatment, species is included
moving_rats = find_rats_that_move(heteros, tags, 8, 9, 3, 4)

# subset tags that never leave a plot
stationary_hets = heteros
plotmovers = unique(moving_rats$tag)
for (i in 1:length(plotmovers)) {
  stationary_hets = subset(stationary_hets, tag != plotmovers[i])
}

# list the stakes it inhabits
tags = unique(stationary_hets$tag)
stakemoves = examine_stake_moves(stationary_hets, tags, 5, 6, 7, 8, 9)

# calculate the distances between each trapping location
plot_stake_moves(stationary_hets, tags, 5, 4, 8, 9)


# plot density of all consecutive movement from rodents within a species 2000-2009
plot(density(pf_meters), main = paste("P. flavus (", length(pftags), " = i, ", length(pf_meters), " = N)", sep = ''), xlab = 'meters', lwd = 2, xlim = c(0,500), ylim = c(0,.1))
plot(density(pp_meters), main = paste("C. penicillatus (", length(pptags), " = i, ", length(pp_meters), " = N)", sep = ''), xlab = 'meters', lwd = 2, xlim = c(0,500), ylim = c(0,.1))
plot(density(pb_meters), main = paste("C. baileyi (", length(pbtags), " = i, ", length(pb_meters), " = N)", sep = ''), xlab = 'meters', lwd = 2, xlim = c(0,500), ylim = c(0,.1))
plot(density(do_meters), main = paste("D. ordii (", length(dotags), " = i, ", length(do_meters), " = N)", sep = ''), lwd = 2, xlim = c(0,500), ylim = c(0,.1))
plot(density(dm_meters), main = paste("D. merriami (", length(dmtags), "= i, ", length(dm_meters), " = N)", sep = ''), xlab = 'meters', lwd = 2, xlim = c(0,500), ylim = c(0,.1))

plot(density(pe_meters), main = paste("P. eremicus (", length(petags), " = i, ", length(pe_meters), " = N)", sep = ''), xlab = 'meters', lwd = 2, xlim = c(0,500), ylim = c(0,.1))
plot(density(pm_meters), main = paste("P. maniculatus (", length(pmtags), " = i, ", length(pm_meters), " = N)", sep = ''), xlab = 'meters', lwd = 2, xlim = c(0,500), ylim = c(0,.1))
plot(density(rm_meters), main = paste("R. megalotis (", length(rmtags), "= i, ", length(rm_meters), " = N)", sep = ''), xlab = 'meters', lwd = 2, xlim = c(0,500), ylim = c(0,.1))

plot(density(sh_meters), main = paste("S. hispidus (", length(shtags), " = i, ", length(sh_meters), " = N)", sep = ''), xlab = 'meters', lwd = 2, xlim = c(0,500), ylim = c(0,.1))
plot(density(sf_meters), main = paste("S. fulviventer (", length(sftags), " = i, ", length(sf_meters), " = N)", sep = ''), xlab = 'meters', lwd = 2, xlim = c(0,500), ylim = c(0,.1))
plot(density(nao_meters), main = paste("N. albigula (", length(naotags), " = i, ", length(nao_meters), " = N)", sep = ''), lwd = 2, xlim = c(0,500), ylim = c(0,.1))

plot(density(ot_meters), main = paste("O. torridus (", length(ottags), " = i, ", length(ot_meters), " = N)", sep = ''), xlab = 'meters', lwd = 2, xlim = c(0,500), ylim = c(0,.1))
plot(density(ol_meters), main = paste("O. leucogaster (", length(oltags), " = i, ", length(ol_meters), " = N)", sep = ''), xlab = 'meters', lwd = 2, xlim = c(0,500), ylim = c(0,.1))


#plot all together
plot(density(pb_meters), main = 'Portal rodent movement', xlab = 'meters', lwd = 2, col = 'peru')
lines(density(pp_meters), col = 'black', lwd = 2)
lines(density(do_meters), col = 'mediumpurple4', lwd = 4, lty = 3)
lines(density(dm_meters), col = 'maroon4', lwd = 4, lty = 3)
lines(density(pf_meters), col = 'hotpink', lwd = 2)

lines(density(pe_meters), col = 'deepskyblue3', lwd = 3, lty = 6)
lines(density(pm_meters), col = 'royalblue4', lwd = 3, lty = 6)
lines(density(rm_meters), col = 'cadetblue', lwd = 3, lty = 6)

lines(density(sh_meters), col = 'indianred', lwd = 3)
lines(density(sf_meters), col = 'brown', lwd = 3)
lines(density(nao_meters), col = 'gray60', lwd = 3)

lines(density(ot_meters), col = 'darkgreen', lwd = 2)
lines(density(ol_meters), col = 'darkolivegreen', lwd = 2)

legend('topright', c('PB','PP','DO','DM','PE','PM','OT','OL'), bty = 'n', lty = c(1,1,3,3,6,6,1,1), lwd = 10, seg.len = 2, 
       col = c('peru', 'black', 'mediumpurple4','maroon4','deepskyblue3', 'royalblue3','darkgreen', 'darkolivegreen'))
abline(v = 70.71, lty = 2, col = 'gray60', lwd = 2)
