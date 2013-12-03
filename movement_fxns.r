# Module containing functions for looking at individual rodent data pertaining to movement


id_unknowns = function(dat, tag_col){
  # give unique numbers to blank tags
  # note: these are 7 digit numbers, so they are longer than any other tag type
  # note: in the Portal data, column 12 is tag, so we are looking for blank or "0" tags to rename
  
unk = 1000000
for (irow in 1:nrow(dat)){
    tag = dat[irow,tag_col]
    unk = unk + 1
    if (tag == "") {
      dat[irow,tag_col] = unk
    }
    else if (tag == "0") {
      dat[irow,tag_col] = unk
    }}
return(dat)
}


find_bad_data = function(dat, tags, sex_col, spp_col){
  # check for consistent sex and species, outputs flagged tags to check
  
flagged_rats = data.frame("tag"=1, "reason"=1, "occurrences"=1)
outcount = 0

for (t in 1:length(tags)){
  tmp <- which(dat$tag == tags[t])
    if (nchar(tags[t]) < 6) {    # is it a pit tag or was it given an unknown tag
      outcount = outcount + 1
    flagged_rats[outcount,] <- c(tags[t], 'not_a_PIT', nrow(dat[tmp,]))
  }
    if (nrow(dat[tmp,]) > 1) {    # if indiv was captured multiple times
    sex_list = dat[tmp,sex_col]
    sex = sex_list[1]
    for (i in 2:length(sex_list)){  # check for consistent sex
      if (sex_list[i] != sex) {
      outcount = outcount + 1
      flagged_rats[outcount,] <- c(tags[t], "sex", nrow(dat[tmp,]))
      break
      }}
    spp_list = dat[tmp,spp_col]
    spp = spp_list[1]
    for (s in 2:length(spp_list)){  # check for consistent species
      if (spp_list[s] != spp){
      outcount = outcount + 1
      flagged_rats[outcount,] <- c(tags[t], "spp", nrow(dat[tmp,]))
      break
      }}
    }}
return(flagged_rats)
}

is_duplicate_tag = function(dat, tags, sex_col, spp_col, tag_col){
  # check the min to max year for a given tag. 
  # If > 4, considered suspicious
  # If multiple species, considered suspicious
  # If adequately resolved, given a new unique tag number, that ends with d for "duplicate"
  numcount = 100
  
  for (t in 1:length(tags)){
    
    if (nchar(tags[t]) < 6){ #only run on ear and toe tags, pit tags are very unlikely to be duplicated
      tmp <- which(dat$tag == tags[t])
      
    if (nrow(dat[tmp,]) > 1) {    # if indiv was captured multiple times
      
      if (max(dat[tmp,1]) - min(dat[tmp,1]) > 4){ #more than 4 years between recaptures?
        
        # check num species recorded? does it look ok if separated on species?
        spp_list = unique(dat[tmp,spp_col])
        for (sp in 1:length(spp_list)) {
          tmp2 = which(dat$tag == tags[t] & dat$species == spp_list[sp])
          if(max(dat[tmp2,1]) - min(dat[tmp2,1]) < 4) {
            newtag = paste(tags[t], numcount, "d", sep = "") #make a new tag to keep separate
            dat[tmp2,tag_col] = newtag
            numcount = numcount + 1
          }
#          else {
#            outcount = outcount + 1
#            flagged_rats[outcount,] <- c(tags[t], "year", nrow(dat[tmp,]))
          }
        }
      }
#    else { #if less than 4 years between captures, check sex
#      
#    }
    }}
  return (dat)
}



find_bad_data2 = function(dat, tags, sex_col, spp_col){
  # check for consistent sex and species, outputs flagged tags to check
  # keeps ear and toe tags, checks for duplicates that appear over too long of a time period (exceed rodent lifetime)
  # or duplicates that appear too close together (multiple indivs with same tag)
  
  flagged_rats = data.frame("tag"=1, "reason"=1, "occurrences"=1)
  outcount = 0
  
  for (t in 1:length(tags)){
    tmp <- which(dat$tag == tags[t])

    if (nrow(dat[tmp,]) > 1) {    # if indiv was captured multiple times
      sex_list = dat[tmp,sex_col]
      sex = sex_list[1]
      for (i in 2:length(sex_list)){  # check for consistent sex
        if (sex_list[i] != sex) {
         outcount = outcount + 1
          flagged_rats[outcount,] <- c(tags[t], "sex", nrow(dat[tmp,]))
          break
        }}
      spp_list = dat[tmp,spp_col]
      spp = spp_list[1]
      for (s in 2:length(spp_list)){  # check for consistent species
        if (spp_list[s] != spp){
          outcount = outcount + 1
          flagged_rats[outcount,] <- c(tags[t], "spp", nrow(dat[tmp,]))
          break
        }}
    }}
  return(flagged_rats)
}
  

subsetDat = function(dataset){
  ## function to subset out proper data 
  ##### will find bad data, then delete it from the dataset and get rid of incompletely sampled periods
  tags = as.character(unique(dataset$tag)) # get list of unique tags
  flags = find_bad_data2(dataset, tags, 10, 9)   # list of flagged data
  #first, mark all uncertain or unmarked sex as "U" for unknown
  badsextags = unique(flags[which(flags$reason == "sex"),1])
  dataset[which(dataset$tag %in% badsextags),10] = "U"
  dataset[which(dataset$sex %in% c("", "P", "Z")),10] = "U" #get rid of other weird typos in sex column
  #get rid of results where we don't know the species for sure
  badspptags = unique(flags[which(flags$reason == "spp"), 1])    
  dataset = dataset[-which(dataset$tag %in% badspptags),] #delete rows where species is unsure
  
  #don't use negative period numbers and periods with only one day of trapping
  dataset = subset(dataset, period != 267 & period != 277 & period != 278 & period != 283 &
                          period != 284 & period != 300 & period != 311 & period != 313 &
                          period != 314 & period != 318 & period != 321 & period != 323 &
                          period != 337 & period != 339 & period != 344 & period != 351)
  return (dataset)
  }



distance_moved = function (data, tags) {
  # Calculatemoved from one time frame to the next, for an individual
  distance = as.numeric()
  # for each individual
  for (t in 1:length(tags)){
    ind_data = data[which(data$tag == tags[t]),] #get data for indiv with tag t
    ind_data = ind_data[order(ind_data$period),] #order chronologically
    # if it was captured more than once
    if (nrow(ind_data) > 1) { 
      for (i in 1:nrow(ind_data)){
        if (i+1 <= nrow(ind_data)){
          meters = sqrt((ind_data[i,8]-ind_data[i+1,8])**2 + (ind_data[i,7]-ind_data[i+1,7])**2)
          distance=append(distance,meters)
        }
      }
    }
  }
  return(distance)
}


sd_avg_mass = function (ind_dat, spp_mean, spp_sd) {
  # Then takes the average mass for an individual (across all recaptures) 
  # and calculates the proportional difference away from that avg. species level mass.

  ind_mean = mean(ind_dat$wgt, na.rm = TRUE)
  
  if (is.na(ind_mean)) {
    # if no mass is recorded (returns NA), make ind_sd = 0 so we don't lose data. 
    #This shouldn't skew the results, but it ISN'T a true zero. Denote appropriately in methods.
    ind_sd = 0
  }
  
  else {
  # num of standard deviations the indiv's mean wgt is from spp mean wgt
  ind_sd = round((ind_mean - spp_mean) / spp_sd , 4) 
  }
  
  return (ind_sd) #number of standard deviations individual is away from capture weight mean
}


feeding_guild = function(speciesname) {
  # grab the species name and decide what feeding guild it is in, based on the lit
  
  # heteromyidae granivores == 1
  if (speciesname %in% list("DO", "DM", "DS", "PB", "PP", "PF", "PH", "PI")) {guild = 1}
  # cricetidae granivores == 2
  else if (speciesname %in% list("PE", "PM", "PL", "RM", "RF", "RO", "BA")) {guild = 2}
  # folivores == 3
  else if (speciesname %in% list("SH", "SF", "SO", "NAO")) {guild = 3 }  
  # insectivores == 4
  else {guild = 4}
  
  return(guild)
}


enumerate_species = function(speciesname) {
  #each species needs to be replaced with a number instead of a name,
  #for input into Program MARK later
  
  if (speciesname == "DO") {speciesnum = 1}
  else if (speciesname == "DM") {speciesnum = 2}
  else if (speciesname == "DS") {speciesnum = 3}
  else if (speciesname == "PB") {speciesnum = 4}
  else if (speciesname == "PP") {speciesnum = 5}
  else if (speciesname == "PF") {speciesnum = 6}
  else if (speciesname == "PE") {speciesnum = 7}
  else if (speciesname == "PM") {speciesnum = 8}
  else if (speciesname == "PL") {speciesnum = 9}
  else if (speciesname == "PH") {speciesnum = 10}
  else if (speciesname == "PI") {speciesnum = 11}
  else if (speciesname == "RM") {speciesnum = 12}
  else if (speciesname == "RF") {speciesnum = 13}
  else if (speciesname == "RO") {speciesnum = 14}
  else if (speciesname == "BA") {speciesnum = 15}
  else if (speciesname == "SH") {speciesnum = 16}
  else if (speciesname == "SF") {speciesnum = 17}
  else if (speciesname == "SO") {speciesnum = 18}
  else if (speciesname == "NAO") {speciesnum = 19}
  else if (speciesname == "OT") {speciesnum = 20}
  else if (speciesname == "OL") {speciesnum = 21}
  
  return(speciesnum)
}

temporal_status = function (speciesname){
  #assigns temporal status to each species based on its across year and within year presence thru span of the entire 30+ year dataset
  #core = 1, intermediate = 2, transient = 3 
  #FIXME: For now, coding as Core vs. Not Core - decide if this is correct later
  
  if (speciesname %in% list("DO", "DM", "PB", "PP", "OT")) {status = 1}
  else if (speciesname %in% list("PE", "RM", "NAO")) {status = 2}
  else if (speciesname %in% list("PF", "PM", "SH", "SF", "OL")) {status = 2}
  
  return(status)
}

noplacelikehome = function (dat, prd, exclosures, breakpoint){
  ### Create a set of MARK capture histories by home vs away from home
  # Creates a movement history to be used in Mark. Matrix is filled in with zeroes (not captured) and later filled in 
  ## 1 (stayed home), and 2 (away from home). 
  ## Home is determined using the mean + 1 sd of the logged data.
  # guild (1 = heteromyid granivore, 2 = cricetid granivore, 3 = folivore, 4 = carnivore)
  # species (assigned using function enumerate_species)
  # status (1 = core, 2 = intermediate, 3 = transient)
  
  tags = unique(dat$tag)
  capture_history = matrix(0, nrow = length(tags), ncol = length(prd))
  covariates = matrix(0, nrow = length(tags), ncol = 7)
    colnames(covariates) = c("male", "female", "unidsex", "sd_mass", "guild", "species", "status")
    group = c(1,2,3) #represent the "group" (core, transient, intermediate)

  # finds the average of all the masses for all non-pregnant adults in the species, as a baseline. 
  # remove juveniles and pregnant females for adult mass estimation
  adult_dat = dat[which(dat$reprod != "J" & dat$pregnant != "P"),] 
  mean_mass = mean(adult_dat$wgt, na.rm = TRUE) 
  sd_mass = sd(adult_dat$wgt, na.rm = TRUE)
  
  # since data is imported by species, we only need to check the first row of data to grab the species name and decide what guild it is in
  # record guild in col 5 of covariates
  covariates[,5] = feeding_guild(dat[1,]$species)
 
  # record species in col 6 of covariates    
  covariates[,6] = enumerate_species(dat[1,]$species)
 
  # record hypothesized status in col 7 of covariates   
  covariates[,7] = temporal_status(dat[1,]$species)

  #loop through each tag to get individual-level data
  for (t in 1:length(tags)) {
    ind_dat = dat[which(dat$tag == tags[t]),] #get data for indiv with tag t
    ind_dat = ind_dat[order(ind_dat$period),] #order chronologically
    
    p1 = min(ind_dat$period) # first capture period
    index = match(p1, prd) # match the period with the index number for the list of periods (will correspond to col num in matrix)
    capture_history[t,index] = 1  #mark first capture with 1 ("home")      
    
    #mark sex in Male, Female, or Unidentified columsn of Covariates
    sex = ind_dat[1,]$sex 
    if (sex == "M") {covariates[t,1] = 1 } 
    else if (sex == "F") {covariates[t,2] = 1 }
    else {covariates[t,3] = 1 }
    
    # record standard deviations away from species average mass as another covariate
    covariates[t,4] = sd_avg_mass(ind_dat, mean_mass, sd_mass) 
    
    for (i in 1:nrow(ind_dat)){ #record capture history data
      
      if (i+1 <= nrow(ind_dat)){
        meters = sqrt((ind_dat[i,8]-ind_dat[i+1,8])**2 + (ind_dat[i,7]-ind_dat[i+1,7])**2)
        pnext = ind_dat[i+1,]$period #next capture period, where the distance will be recorded in the matrix
        
        if (meters <= breakpoint) {
          dist = 1 #captured close to "home"
        }
        
        else if (meters > breakpoint) {
          dist = 2 #captured far from "home"
        }

        if (ind_dat[i+1,]$plot %in% exclosures){ #was it captured on an exclosure?
          covariates[t,group] = covariates[t,group]*-1
        }
        
        index = match(pnext, prd)
        capture_history[t,index] = dist #mark subsequent captures 
      }
    }  
  }
  mark_df = concat_ch(capture_history, covariates)
  return(mark_df)
}


concat_ch = function (ch_matrix, cov_matrix){
  #concatenates columns representing capture histories, to be used in later MARK analyses
  ch_df <- data.frame(ch_matrix)
  encounters <- do.call(paste, c(ch_df[c(names(ch_df))], sep = '')) # makes a vector of all the capture histories
      semicol <- rep(";", nrow(ch_df)) # makes a vector of semicolons
      mark_df <- cbind(encounters, cov_matrix, semicol) # binds the capture, cov, and semicolon data together into a dataframe
  return (mark_df)
}


mean_win_yr_occ = function (data){
  #finds the mean within year occupancy for each month for a given species, returns a single value
  uniq_mos = c(10, 8, 10, 10, 8, 9, 9, 10, 12, 12) #not all months were trapped in all years
  mos = c(1:12)
  years = c(2000:2009)
  
  proportion_mos = c()
  
  for (y in 1:length (years)){
    yr_data = subset(data, yr == years[y])
    if(nrow(yr_data) > 0) {  #don't use years where it wasn't captured
      m = length(unique(yr_data$mo))/uniq_mos[y]
      proportion_mos = append(proportion_mos, m)
    }
  }
    months = round(mean(proportion_mos),4)
    
  return (months)
}


mean_win_yr_alldat = function (data){
  #TODO: uniq_mos will not be correct now that we have changed the number of years sampled
    # add a fxn that will concatenate the number of unique mos sampled in a given year
  #finds the mean within year occupancy for each month for a given species, returns a single value
  uniq_mos = c(12, 10, 12, 11, 10, 12,  9, 11, 10, 11, 12, 11, 12,
               12, 10, 12, 11, 10, 12, 10, 11, 10, 10, 11, 11, 10,
               12, 11, 11, 12, 12, 12, 8, 9, 12) #not all months were trapped in all years
  mos = c(1:12)
  years = c(1978:2012)
  
  proportion_mos = c()
  
  for (y in 1:length (years)){
    yr_data = subset(data, yr == years[y])
    if(nrow(yr_data) > 0) {  #don't use years where it wasn't captured
      m = length(unique(yr_data$mo))/uniq_mos[y]
      proportion_mos = append(proportion_mos, m)
    }
  }
  months = round(mean(proportion_mos),4)
  
  return (months)
}


mean_mo_repro = function (femaledata){
  #returns the proportion of females that are reproductive (nipples enlarged - E - red - R- or both - B - or pregnant - P-)
  #on average in a given month across all the years. Only looks at data during years and months in which the species is present. 
  mos = c(1:12)
  years = sort(unique(femaledata$yr)) #only look at data during years in which the species is present

  avg_r_mo = c()
  
  for (m in 1:length(mos)){
    mo_repros = c()
    
    for (y in 1:length(years)){
      tmp = subset(femaledata, yr == years[y])
      tmp = subset(tmp, mo == mos[m])
      
      if (nrow(tmp) > 0){
        num_females = nrow(tmp)
        repro = subset(tmp, nipples == "E" | nipples == "B" | nipples == "R" | pregnant == "P")
        prop_repro = nrow(repro)/num_females
        mo_repros = append(mo_repros, prop_repro)
      }
      else {
        mo_repros = append(mo_repros, NA)
      }
    }
    avg_mo = round(mean(mo_repros, na.rm = T),4)
    avg_r_mo = append(avg_r_mo, avg_mo)
  }
  
  return(avg_r_mo)
}


mo_repro = function (femaledata){
  #returns the proportion of females that are reproductive (nipples enlarged - E - red - R- or both - B - or pregnant - P-)
  #in each year. Only looks at data during years and months in which the species is present. 
  mos = c(1:12)
  years = sort(unique(femaledata$yr)) #only look at data during years in which the species is present
  species = femaledata[1,9]
  
  r_mo_df =data.frame("year" = 1, "month" = 1, "proprepro" = 1, "numfemales" = 1, "species" = 1)
  
  for (m in 1:length(mos)){
      tmp = subset(femaledata, mo == mos[m])
      
    for (y in 1:length(years)){
      tmp = subset(tmp, yr == years[y])
      
      if (nrow(tmp) > 0){
        num_females = nrow(tmp)
        repro = subset(tmp, nipples == "E" | nipples == "B" | nipples == "R" | pregnant == "P")
        prop_repro = round(nrow(repro)/num_females, 4)
        mo_repros = c(years[y], mos[m], prop_repro, num_females, species)
      }
      else {
        num_females = 0
        mo_repros = c(years[y], mos[m], NA, num_females, species)
      }
      r_mo_df = rbind(r_mo_df, mo_repros)
    }
  }
  
  return(r_mo_df[-1,])
}


count_repro = function(reprodata){
  #input a df of reproductive data. Outputs the number of times that the individual was "uniquely" reproductive.
  # unique reproduction is defined by having a break between reproductive status in the trapping record.
 
  rh = 1   #reproductive history counter. Starts at one. Increment when the difference between reproductive events is > 1
  
  if (nrow(reprodata) > 1) { 
    prds = reprodata$period
    
    for (p in 1:length(prds)){
      
      if ((p + 1) <= length(prds)){ #can't go past the number of periods that exist for the individual
        diff = prds[p+1] - prds[p] #subtract the periods
        
        if (diff > 1){
          rh = rh + 1 #increment by one if the difference between reproductive events is > 1 period
        }
      }
    }  
  }
  return (rh)
}

indiv_repro = function (femaledata){
  # returns the reproductive history of females in the time series
  # for each individual, all rows of data where nipples are red, enlarged or both, or it is pregnant,
  # are set aside. Within a year, the trapping periods are ordered to detect timing of reproduction
  # and the number of unique reproduction events, per year, for each individual female.
  tags = unique(femaledata$tag)

  reprod_df = data.frame("tag" = 1, "year" = 1, "num_reprod" = 1)
  
  for (t in 1:length(tags)) {
    indiv = subset(femaledata, tag == tags[t]) #get individual data, sorted by chronologically by period
    indiv = indiv[order(indiv$period),] #order chronologically
    
    years = sort(unique(indiv$yr))
    
    for (y in 1:length(years)){
      tmp = subset(indiv, yr == years[y])
      repro = subset(tmp, nipples == "E" | nipples == "B" | nipples == "R" | pregnant == "P")
      
      if (nrow(repro) > 0){
        numreprod = count_repro(repro)
      }
      
      else {
        numreprod = 0
      }
      
    }
    data = c(tags[t], years[y], numreprod)
    reprod_df = rbind(reprod_df, data)
  }
  return (reprod_df[-1,])
}


allyrs_abun = function(sp_data, years){  #TODO: fix inputs so that you can input a range of years to use, instead of assuming fixed range
  #function to find abundance in years in which the species is present. input the species dataframe, For each
  # year the species occurs in, calculates the abundance (total number of unique individuals). Returns a vector
  years = years #TODO: Make sure inputs in rodent_wrapper.r are correct
  abun = c()
  
  for (y in 1:length(years)){
    dat = subset(sp_data, sp_data[,1] == years[y])
    if (length(dat) > 0) {
      indivs = sort(unique(dat$tag))
      abun = append(abun, length(indivs))
    }
  }
  return (abun)
}
