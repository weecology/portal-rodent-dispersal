# for manipulating  rodent movement data in MARK
## analyze survival and dispseral probabilities for species and guilds

library(RMark) #can't install on mac? Update R version?

#mac access
wd = "/Users/sarah/Desktop/portal-rodent-dispersal/"
setwd(wd)

# bring in the inp file and convert it to RMark format (.txt needs to be .inp, change manually)
MSdata <- convert.inp("do_mark.inp",
                      group.df=data.frame(sex=c("male","female","unidsex")))

# Build up the model. Looking at sex effects on dispersal/survival
MSprocess <- process.data(MSdata,model="Multistrata",begin.time=2000,
                           group=c("sex"))

MS.ddl <- make.design.data(MSprocess)


