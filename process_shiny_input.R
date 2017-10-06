library('reshape2')
library('plyr')
library('dplyr')

#set the working directory using 'setwd()'

filename <- "Simulated1"
secperframe <- 0.1
pixelsize <- 71 #nm

msd <- read.csv(paste0("Output/",filename,"_msdfull_R.csv"), header = FALSE)
ids <- read.csv(paste0("Output/",filename,"_msdids_R.csv"), header = TRUE)
fits <- read.csv(paste0("Output/",filename,"_fit10_inter22__new_R.csv"), header = TRUE)
pos <- read.csv(paste0("Output/",filename,"_trajpos_R.csv"), header = TRUE)

pixelsizeum <- pixelsize/1000

msd.m <- melt(msd, na.rm = TRUE, id.vars = 'V1')
msd.m.rn <- msd.m %>% mutate(rownum = as.integer(substring(as.character(variable),2)))

data1 <- left_join(msd.m.rn,ids,by = 'rownum') %>% select(ID = id, SizeClass = sizeclass, lagframe = V1, MSD = value, length) 
data1$SizeClass <- factor(data1$SizeClass)
data1$SizeClass <- revalue(data1$SizeClass,c("2"="L", "1"="S"))
data1$SizeClass <- as.character(data1$SizeClass)

fits$SizeClass <- as.character(fits$SizeClass)

data <- left_join(data1,fits, by = c("ID","SizeClass")) %>% mutate(uniqueID = paste0(as.character(ID),SizeClass)) %>%
  mutate(lagtime = lagframe*secperframe, MSDum = MSD*pixelsizeum^2, lengthtime = length*secperframe) %>% select(-Length,-Size,-NumSizeErrors)

positions = pos
positions$class <-  revalue(as.character(positions$class),c("2"="L", "1"="S"))
positions <- positions %>% mutate(uniqueID = paste0(as.character(id),class), realtime = time * secperframe) %>% 
  select(realtime,frametime = time,x,y,uniqueID,ID = id, SizeClass = class)

#add MSD
rm(list=setdiff(ls(), c("data","positions","filename","pixelsize","secperframe")))
v <- 0.2
save.image(file = paste0("Output/",filename,"_shiny.Rdata"))
