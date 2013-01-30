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
    if (nchar(tags[t]) < 6) {
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


subsetDat = function(dataset){
  ## function to subset out proper data 
  ##### will find bad data, then delete it from the dataset and get rid of incompletely sampled periods
  tags = as.character(unique(dataset$tag)) # get list of unique tags
  flags = find_bad_data(dataset, tags, 10, 9)   # list of flagged data
  badtags = unique(flags$tag)    # get list of unique 'bad tags'
  for (i in 1:length(badtags)) {  #deletes bad tags from dataset
    dataset = subset(dataset, tag != badtags[i])
  }  
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


sd_avg_mass = function (dat, ind_dat) {
  # finds the average of all the masses for all the captures in the species data, as a baseline. 
  ## Then takes the average mass for an individual and calculates the proportional difference away from that avg. species
  ## level mass.
  spp_mean = mean(dat$wgt, na.rm = TRUE) #include ALL indivs (should I NOT INCLUDE juveniles and pregnant indivs?)
  spp_sd = sd(dat$wgt, na.rm = TRUE)
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


noplacelikehome = function (dat, prd, exclosures, breakpoint){
  ### Create a set of MARK capture histories by home vs away from home
  # Creates a movement history to be used in Mark. Matrix is filled in with zeroes (not captured) and later filled in 
  ## 1 (stayed home), and 2 (away from home). 
  ## Home is determined using the mean + 1 sd of the data.
  
  tags = unique(dat$tag)
  capture_history = matrix(0, nrow = length(tags), ncol = length(prd))
  covariates = matrix(0, nrow = length(tags), ncol = 7)
    colnames(covariates) = c("male", "female", "unidsex", "sd_mass", "hgran", "cgran", "foli")
    group = c(1,2,3) #represent the "group"
  
  # record guild in dummy variables (cols 5:7)
  # since data is imported by species, we only need to check the first row of data to grab the species name and decide what guild it is in
  if (dat[1,]$species %in% list("DO", "DM", "PB", "PP", "PF")){ 
    covariates[, 5] = 1 }
  else if (dat[1,]$species %in% list("PE", "PM", "RM")){
    covariates[, 6] = 1}
  else if (dat[1,]$species %in% list("SH", "SF", "NAO")){
    covariates[,7] = 1 }  #all zeros indicate insectivores    
  
  #loop through each tag to get individual-level data
  for (t in 1:length(tags)) {
    ind_dat = dat[which(dat$tag == tags[t]),] #get data for indiv with tag t
    ind_dat = ind_dat[order(ind_dat$period),] #order chronologically
    
    p1 = min(ind_dat$period) # first capture period
    index = match(p1, prd) # match the period with the index number for the list of periods (will correspond to col num in matrix)
    capture_history[t,index] = 1  #mark first capture with 1 ("home")      
    
    sex = ind_dat[1,]$sex # don't need to acct for disputes in sex b/c should be already deleted (in flagged data fxn)
    if (sex == "M") {
      covariates[t,1] = 1 } 
    else if (sex == "F") {
      covariates[t,2] = 1 }
    else { 
      covariates[t,3] = 1 }
    
    sd_mass = sd_avg_mass(dat, ind_dat) # record standard deviations away from species average mass
    covariates[t,4] = sd_mass 
    
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
  mos = c(1:12)
  years = sort(unique(data$yr)) #only count during years in which the species is present
  
  proportion_mos = c()
  
  for (y in 1:length (years)){
    yr_data = subset(data, yr == years[y])
    m = length(unique(yr_data$mo))/12
    proportion_mos = append(proportion_mos, m)
  }
  months = round(mean(proportion_mos),4)
  
  return (months)
}

