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
  # num of standard deviations the indiv's mean wgt is from spp mean wgt
  ind_sd = round((ind_mean - spp_mean) / spp_sd , 4) 
  
  return (ind_sd) #number of standard deviations individual is away from capture weight mean
}


noplacelikehome = function (dat, prd, exclosures, breakpoint){
  ### Create a set of MARK capture histories by home vs away from home
  # Creates a movement history to be used in Mark. Matrix is filled in with zeroes (not captured) and later filled in 
  ## 1 (stayed home), and 2 (away from home). 
  ## Home is determined using the mean + 1 sd of the data.
  
  tags = unique(dat$tag)
  MARK_distance = matrix(0, nrow = length(tags), ncol = length(prd) + 4)
    group = c(ncol(MARK_distance)-3, ncol(MARK_distance)-2, ncol(MARK_distance)-1) #represent the "group"
  
  for (t in 1:length(tags)) {
    ind_dat = dat[which(dat$tag == tags[t]),] #get data for indiv with tag t
    ind_dat = ind_dat[order(ind_dat$period),] #order chronologically
    
    p1 = min(ind_dat$period) # first capture period
    index = match(p1, prd) # match the period with the index number for the list of periods (will correspond to col num in matrix)
    MARK_distance[t,index] = 1  #mark first capture with 1 ("home")      
    
    sex = ind_dat[1,]$sex # don't need to acct for disputes in sex b/c should be already deleted (in flagged data fxn)
    if (sex == "M") {
      MARK_distance[t,ncol(MARK_distance)-3] = 1 } 
    else if (sex == "F") {
      MARK_distance[t,ncol(MARK_distance)-2] = 1 }
    else { 
      MARK_distance[t,ncol(MARK_distance)-1] = 1 }
    
    sd_mass = sd_avg_mass(dat, ind_dat)
    MARK_distance[t,ncol(MARK_distance)] = sd_mass #put sd of mass in last col
    
    for (i in 1:nrow(ind_dat)){
      
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
            MARK_distance[t,group] = MARK_distance[t,group]*-1
        }
        
        index = match(pnext, prd)
        MARK_distance[t,index] = dist #mark subsequent captures 
      }
    }  
  }
  mark_df = concat_ch(MARK_distance)  
  return(mark_df)
}

concat_ch = function (matrix){
  #concatenates columns representing capture histories, to be used in later MARK analyses
  # A is a placeholder for the concatenated data (code from D. Koons)
  A <- data.frame(matrix)
  concat <- paste(A$X1,A$X2,A$X3,A$X4,A$X5,A$X6,A$X7,A$X8,A$X9,A$X10,A$X11,A$X12,
                   A$X13,A$X14,A$X15,A$X16,A$X17,A$X18,A$X19,A$X20,A$X21,A$X22,A$X23,A$X24,A$X25,
                   A$X26,A$X27,A$X28,A$X29,A$X30,A$X31,A$X32,A$X33,A$X34,A$X35,A$X36,A$X37,A$X38,A$X39,
                   A$X40,A$X41,A$X42,A$X43,A$X44,A$X45,A$X46,A$X47,A$X48,A$X49,A$X50,A$X51,A$X52,
                   A$X53,A$X54,A$X55,A$X56,A$X57,A$X58,A$X59,A$X60,A$X61,A$X62,A$X63,A$X64,A$X65,
                   A$X66,A$X67,A$X68,A$X69,A$X70,A$X71,A$X72,A$X73,A$X74,A$X75,A$X76,A$X77,A$X78,
                   A$X79,A$X80,A$X81,A$X82,A$X83,A$X84,A$X85,A$X86,A$X87,A$X88,A$X89,A$X90,A$X91,
                   A$X92,A$X93,A$X94,A$X95,A$X96,A$X97,A$X98,A$X99,A$X100,A$X101,A$X102,A$X103,A$X104,
                   A$X105,A$X106,A$X107,A$X108,A$X109,A$X110,A$X111,A$X112,A$X113,A$X114,A$X115,
                   A$X116,A$X117,A$X118,A$X119,A$X120,sep='')
      semicol <- rep(";", nrow(A))
      mark_df <- cbind(concat,A$X121,A$X122,A$X123,A$X124,semicol)
  return (mark_df)
}

format_wrt_mark = function (dataframe, pathname){
  # format and write files, includes concatenating recapture histories
}




#################### OLD FUNCTIONS ? TO FIX AND WORK ON...

### Create a set of capture histories by treatment type
create_trmt_hist = function(dat, tags, prd){

MARK_data = data.frame("captures"=1, "censored"=1, "tags"=1, "species"=1, "sex"=1, "avg_weight"=1)
outcount = 0

for (t in 1:length(tags)){
  capture_history = "" #create empty string
  for (p in 1:length(prd)){
    tmp<-which(dat$tag==tags[t] & dat$period==prd[p])
    if (nrow(dat[tmp,]) == 0) {
      state = "0"
      capture_history = paste(capture_history, state, sep="")}
    else if (nrow(dat[tmp,] > 0)) {
      if (dat[tmp,17]==1){
        state = "A"
        capture_history = paste(capture_history, state, sep="")}
      else if (dat[tmp,17]==2){
        state = "B"
        capture_history = paste(capture_history, state, sep="")}
      else if (dat[tmp,17]==3){
        state = "C"
        capture_history = paste(capture_history, state, sep="")}
      }
    }
    tmp2<-which(dat$tag==tags[t])
    censored = 1
    for (irow in nrow(dat[tmp2,])){
    if (dat[tmp2,17] == 3) {
      censored = -1
      break}}
    spp = unique(dat[which(dat$tag==tags[t]),6])
    sex = unique(dat[which(dat$tag==tags[t]),7])
    avg_mass = mean(dat[which(dat$tag==tags[t]),11])
    outcount = outcount + 1
    MARK_data[outcount,] <- c(capture_history, censored, tags[t], spp, sex, avg_mass)
  }
return(MARK_data)
}

## Create a set of capture histories by plot
create_plot_hist = function(dat, tags, prd) {
MARK_data_by_plot = data.frame("captures"=1, "censored"=1, "tags"=1, "species"=1, "sex"=1, "avg_weight"=1)
outcount = 0

for (t in 1:length(tags)){
  capture_history = "" #create empty string
  for (p in 1:length(prd)){
    tmp<-which(dat$tag==tags[t] & dat$period==prd[p])
    if (nrow(dat[tmp,]) == 0) {
      state = "0"
      capture_history = paste(capture_history, state, sep=",")}
    else if (nrow(dat[tmp,] > 0)) {
      state = as.character(dat[tmp,4])
      capture_history = paste(capture_history, state, sep=",")}
      }
    tmp2<-which(dat$tag==tags[t])
    censored = 1
    for (irow in nrow(dat[tmp2,])){
    if (dat[tmp2,17] == 3) {
      censored = -1
      break}}
    spp = unique(dat[which(dat$tag==tags[t]),6])
    sex = unique(dat[which(dat$tag==tags[t]),7])
    avg_mass = mean(dat[which(dat$tag==tags[t]),11])
    outcount = outcount + 1
    MARK_data_by_plot[outcount,] <- c(capture_history, censored, tags[t], spp, sex, avg_mass)
  }
return(MARK_data_by_plot)
}

# Find rats that are moving around
find_rats_that_move = function(dat, tags, spp_col, sex_col, trmt_col, plot_col){
moving_rats = data.frame("tag"=1, "spp"=1, "sex"=1, "loc"=1, "occurrences"=1)
outcount = 0

for (t in 1:length(tags)){
  tmp <- which(dat$tag == tags[t])
  spp = unique(dat[tmp,spp_col])
  sex = unique(dat[tmp,sex_col])
  if (nrow(dat[tmp,]) > 1) {
    trmt_list = dat[tmp,trmt_col]
    trmt = trmt_list[1]
    for (i in 2:length(trmt_list)){
      if (trmt_list[i] != trmt) {
      outcount = outcount + 1
      moving_rats[outcount,] <- c(tags[t], spp, sex, "trmt", nrow(dat[tmp,]))
      break
      }}
    plot_list = dat[tmp,plot_col]
    plot = plot_list[1]
    for (p in 2:length(plot_list)){
      if (plot_list[p] != plot){
      outcount = outcount + 1
      moving_rats[outcount,] <- c(tags[t], spp, sex, "plot", nrow(dat[tmp,]))
      break
      }}
    }}
return(moving_rats)
}

# find num captures/rat
num_captures = function(dat, tags){

rat_catches = data.frame("tag"=1, "spp"=1, "sex"=1, "captures"=1)
outcount = 0

for (t in 1:length(tags)){
  tmp <- which(dat$tag == tags[t])
  spp = unique(dat[tmp,6])
  sex = unique(dat[tmp,7])
  captures = nrow(dat[tmp,])
  outcount = outcount + 1
  rat_catches[outcount,] <- c(tags[t], spp, sex, captures)
  }
rat_catches$captures = as.numeric(rat_catches$captures)
return(rat_catches)
}

### Creates a dataframe to look at movement
examine_trmt_moves = function(dat, tags){

moves_trmt = data.frame("tags"=1, "sp"=1, "sex"=1, "num_moves"=1, "move_list"=1, "censored"=1,
                        "c2r"=1, "r2c"=1, "c2e"=1, "r2e"=1)
outcount = 0

for (t in 1:length(tags)){
  move_list = as.character() #create empty string
  ind_dat=dat[which(dat$tag==tags[t]),]  # subset data for individual with tag t
  if (nrow(ind_dat) > 1) {  # only record data for indivs with multiple captures
    for (i in 1:nrow(ind_dat)){
      if (ind_dat[i,17]==1){
        state = "A"
        move_list <- append(move_list, state)
        censored = 1 }
      else if (ind_dat[i,17]==2){
        state = "B"
        move_list <- append(move_list, state)
        censored = 1 }
      else if (ind_dat[i,17]==3){
        state = "C"
        move_list <- append(move_list, state)
        censored = -1 }
      }
    type = move_list[1]
    moves = 0
    new_list = ""
    for (l in 1:length(move_list)){
      new_list = paste(new_list, move_list[l], sep="")
      c2r = 0
      r2c = 0
      c2e = 0
      r2e = 0
      if (move_list[l]!=type) {
        moves = moves + 1
        if (type == "A" &  move_list[l] == "B") {
          c2r = c2r + 1}
        else if (type == "A" & move_list[l] == "C") {
          c2e = c2e + 1}
        else if (type == "B" & move_list[l] == "A") {
          r2c = r2c + 1}
        else if (type == "B" & move_list[l] == "C") {
          r2e = r2e + 1}
        type = move_list[l] 
        }}
    spp = unique(ind_dat[,6])
    sex = unique(ind_dat[,7])
    outcount = outcount + 1
    moves_trmt[outcount,] <- c(tags[t], spp, sex, moves, new_list, censored, c2r, r2c, c2e, r2e)
  }}
return(moves_trmt)
}

### Creates a dataframe to look at movement among plots
examine_plot_moves = function(dat, tags){

moves_plot = data.frame("tags"=1, "sp"=1, "sex"=1, "num_moves"=1, "neighborhood"=1, "move_list"=1)
outcount = 0

for (t in 1:length(tags)){
  move_list = as.numeric() #create empty string
  ind_dat=dat[which(dat$tag==tags[t]),]  # subset data for individual with tag t
  if (nrow(ind_dat) > 1) {  # only record data for indivs with multiple captures
    for (i in 1:nrow(ind_dat)){
      state = ind_dat[i,4]
      move_list <- append(move_list, state)}
      }
    type = move_list[1]
    neighbor = ""
    moves = 0
    new_list = ""
    for (l in 1:length(move_list)){
      new_list = paste(new_list, move_list[l], sep=",")
      if (move_list[l]!=type) {
        moves = moves + 1
        near = is.neighbor(move_list[l], type) 
        neighbor = paste(neighbor, near, sep = "")
        type = move_list[l] }}
    spp = unique(ind_dat[,6])
    sex = unique(ind_dat[,7])
    outcount = outcount + 1
    moves_plot[outcount,] <- c(tags[t], spp, sex, moves, neighbor, new_list)
  }
return(moves_plot)
}

# determine if a plot is a neighbor
is.neighbor = function(plotA, plotB){

if (plotA == 1) {
  if (plotB == 2 | plotB == 7) { neighbor = "Y" }
  else { neighbor = "N" }}
else if (plotA == 2) {
  if (plotB == 1 | plotB == 3 | plotB == 8 ) { neighbor = "Y" }
  else { neighbor = "N" }}
else if (plotA == 3) {
  if (plotB == 2 | plotB == 4 | plotB == 8 | plotB == 9 | plotB==10) { neighbor = "Y" }
  else { neighbor = "N" }}
else if (plotA == 4) {
  if (plotB == 3 | plotB == 5 | plotB == 9 | plotB == 10 | plotB == 11) { neighbor = "Y" }
  else { neighbor = "N" }}
else if (plotA == 5) {
  if (plotB == 4 | plotB == 6 | plotB == 11 | plotB == 12 | plotB == 24) { neighbor = "Y" }
  else { neighbor = "N" }}
else if (plotA == 6) {
  if (plotB == 5 | plotB==24 | plotB == 12) { neighbor = "Y" }
  else { neighbor = "N" }}
else if (plotA == 7) {
  if (plotB == 1 | plotB == 2 | plotB == 13 | plotB == 8 ) {neighbor = "Y" }
  else { neighbor = "N" }}
else if (plotA == 8 ) {
  if (plotB == 2 | plotB == 7 | plotB == 9 | plotB == 14 ) { neighbor = "Y"}
  else { neighbor = "N" }}
else if (plotA == 9 ) {
  if (plotB == 8 | plotB == 10 | plotB==3 | plotB == 4 | plotB == 16 | plotB == 15) { neighbor = "Y" }
  else { neighbor = "N" }}
else if (plotA == 10 ) {
  if (plotB == 9 | plotB == 11 | plotB == 16 | plotB == 4 | plotB == 17 ) { neighbor = "Y" }
  else { neighbor = "N" }}
else if (plotA == 11 ) {
  if (plotB == 10 | plotB == 12 | plotB == 5 | plotB == 16 | plotB == 17 | plotB == 18) { neighbor = "Y" }
  else { neighbor = "N" }}
else if (plotA == 12 ) {
  if (plotB == 11 | plotB == 5 | plotB == 6 | plotB == 24 | plotB == 19 | plotB == 18) { neighbor = "Y" }
  else { neighbor = "N" }}
else if (plotA == 13) {
  if (plotB == 7 | plotB == 14 | plotB == 8) { neighbor = "Y" }
  else { neighbor = "N" }}
else if (plotA == 14) {
  if (plotB == 13 | plotB == 8 | plotB == 15 | plotB == 7) {neighbor = "Y" }
  else { neighbor = "N" }}
else if (plotA == 15) {
  if (plotB == 14 | plotB == 16 | plotB == 9 | plotB == 8 | plotB == 20) { neighbor = "Y" }
  else { neighbor = "N" }}
else if (plotA == 16) {
  if (plotB == 15 | plotB == 17 | plotB == 9 | plotB == 10 | plotB == 21 | plotB == 20) { neighbor = "Y" }
  else { neighbor = "N" }}
else if (plotA == 17) {
  if (plotB == 16 | plotB == 18 | plotB == 10 | plotB == 11 | plotB == 21 | plotB == 22) { neighbor = "Y" }
  else { neighbor = "N" }}
else if (plotA == 18) {
  if (plotB == 17 | plotB == 19 | plotB == 12 | plotB == 22 | plotB == 23) { neighbor = "Y" }
  else { neighbor = "N" }}
else if (plotA == 19) {
  if (plotB == 12 | plotB == 18 | plotB == 24) { neighbor = "Y" }
  else { neighbor = "N" }}
else if (plotA == 20) {
  if (plotB == 15 | plotB == 16 | plotB == 21) { neighbor = "Y" }
  else { neighbor = "N" }}
else if (plotA == 21) {
  if (plotB == 20 | plotB == 22 | plotB == 16 | plotB == 15) { neighbor = "Y" }
  else { neighbor = "N" }}
else if (plotA == 22) {
  if (plotB == 21 | plotB == 23 | plotB == 17 | plotB == 18) { neighbor = "Y" }
  else { neighbor = "N" }}
else if (plotA == 23) {
  if (plotB == 22 | plotB == 18) { neighbor = "Y" }
  else { neighbor = "N" }}
else if (plotA == 24) {
  if (plotB == 6 | plotB == 19 | plotB == 12) { neighbor == "Y" }
  else { neighbor = "N" }}
else { neighbor = "ERROR" }

return (neighbor)
}

### Creates a dataframe to look at movement among plots          # FIX ME   #WHY ARE THERE MULTIPLE CAPS PER PERIOD?? FIX!
examine_stake_moves = function(dat, tags, stake_col, e_col, n_col, spp_col, sex_col){

moves_stake = data.frame("tags"=1, "sp"=1, "sex"=1, "avg_mass"=1, "num_moves"=1, "stake_list"=1)
outcount = 0
#distances = as.numeric()

for (t in 1:length(tags)){
  stake_list = as.numeric() #create empty string
  ind_dat=dat[which(dat$tag==tags[t]),]  # subset data for individual with tag t
  ind_dat = ind_dat[order(ind_dat[,2]),] # order by period number so time is sequential
  if (nrow(ind_dat) > 1) {  # only record data for indivs with multiple captures
    for (i in 1:nrow(ind_dat)){
      stake = ind_dat[i,stake_col]
      stake_list <- append(stake_list, stake)}  
    type = stake_list[1]
    moves = 0
    new_stake_list = ""
    for (l in 1:length(stake_list)){
      new_stake_list = paste(new_stake_list, stake_list[l], sep=" ")
      if (stake_list[l]!=type) {
        moves = moves + 1
        type = stake_list[l] }}
    if (moves > 1) {  
#    dist_list = as.numeric()
#    for (i in 1:nrow(ind_dat)){
#      if (i > 1){
#      dist = round(UTM_dist(ind_dat[i,e_col], ind_dat[i-1,e_col], ind_dat[i,n_col], ind_dat[i-1,n_col]),3)
#      dist_list = append(dist_list, dist)
#      }}
    avg_mass = averageMass(ind_dat)
    spp = unique(ind_dat[,spp_col])
    sex = unique(ind_dat[,sex_col])
    outcount = outcount + 1
    moves_stake[outcount,] <- c(tags[t], spp, sex, avg_mass, moves, new_stake_list)
  }}}
return(moves_stake)
}

# calculte dist between two UTMs, dist=sqrt((x2-x1)^2 + (y2-y1)^2) 
UTM_dist = function (e1, e2, n1, n2){
  dist = sqrt((e2-e1)^2 + (n2-n1)^2)
  return(dist)
}

# Get average mass from all rows of data pertaining to individual, don't count Juvenile or pregnant mass
averageMass = function(data) {
  mass_data = subset(data, reprod != "J" & pregnant != "P")
  masses = subset(mass_data, mass > 0, na.rm=T)
  avg_mass = sum(masses$mass)/length(masses$mass)
  return(round(avg_mass,3))}
                           
# break data up into x and y coordinates
xy_data = function (stake) {

if (stake==71|stake==72|stake==73|stake==74|stake==75|stake==76|stake==77) {
    x = 7  }
else if (stake==61|stake==62|stake==63|stake==64|stake==65|stake==66|stake==67) {
    x = 6 }
else if (stake==51|stake==52|stake==53|stake==54|stake==55|stake==56|stake==57) {
    x = 5 }
else if (stake==41|stake==42|stake==43|stake==44|stake==45|stake==46|stake==47) {
    x = 4 }
else if (stake==31|stake==32|stake==33|stake==34|stake==35|stake==36|stake==37) {
    x = 3 }
else if (stake==21|stake==22|stake==23|stake==24|stake==25|stake==26|stake==27) {
    x = 2 }
else { x = 1}
if (stake==71|stake==61|stake==51|stake==41|stake==31|stake==21|stake==11) {
    y = 1 }
else if (stake==72|stake==62|stake==52|stake==42|stake==32|stake==22|stake==12) {
    y = 2 }
else if (stake==73|stake==63|stake==53|stake==43|stake==33|stake==32|stake==31) {
    y = 3 }
else if (stake==74|stake==64|stake==54|stake==44|stake==34|stake==24|stake==14) {
    y = 4 }
else if (stake==75|stake==65|stake==55|stake==45|stake==35|stake==25|stake==15) {
    y = 5 }
else if (stake==76|stake==66|stake==56|stake==46|stake==36|stake==26|stake==16) {
    y = 6 }
else  { y = 7 }
xy=c(x,y)
return(xy)
}

# plot stake moves, takes data from those in one plot that are recaptured at least three times
plot_stake_moves = function (dat, tags, stake_col, plot_col, spp_col, sex_col) {
for (t in 1:length(tags)){
  coords = data.frame("x"=1,"y"=1)
  outcount = 0
  ind_dat=dat[which(dat$tag==tags[t]),]  # subset data for individual with tag t
  ind_dat = ind_dat[order(ind_dat[,2]),] # order by period number so time is sequential
  if (nrow(ind_dat) > 2) {  # only record data for indivs with multiple captures
    for (i in 1:nrow(ind_dat)){
      stake = ind_dat[i,stake_col]
      xy = xy_data(stake)
      outcount = outcount + 1
      coords[outcount,]<-xy
      } 
    avg_mass = averageMass(ind_dat)
    spp = unique(ind_dat[,spp_col])
    sex = unique(ind_dat[,sex_col])
    plot_num = unique(ind_dat[,plot_col])
    plot(coords$x, coords$y, xlim = c(1,7), ylim = c(1,7), type = "b", 
          main = paste(avg_mass, "g", sex, spp, "in plot", plot_num, ",recaptures =", nrow(coords), sep = " "))
    points(coords[1,1], coords[1,2], pch = 17, col = "indianred") #start
    points(coords[nrow(coords),1],coords[nrow(coords),2],pch = 15, col = "blue") #stop
  }}}


# plot number of individuals captured in each period (since each record has its own row, we can just count the number
### of occurences of each period)
plot_freq_by_prd = function(data, name) {
  require(plyr)
  
  prdcount = count(data, vars = "period")
  plot(prdcount$period, prdcount$freq, type = "b", pch = 19, xlab = "Period", 
       ylab = paste("num individual", name, "captured", sep = ' '))
}


#plot number of recaptures as a histogram for each species
plot_recap_hist = function (data, name) {
  recapcount = count(data, vars = "tag")
  recaps = recapcount$freq
  bw = seq(min(recaps), max(recaps), 1)
  hist(recaps, breaks = bw, xlab = paste("num", name, "recaptures", sep = " "), ylab = "count",
       col = "gray60", border = "white", main = "")
}


