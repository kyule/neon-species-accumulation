### This is the primary file for analyzing NEON species accumulation data
### Want to understand years of sampling required, but sampling differs across years, need to scale by the number of traps


#### BEFORE RUNNING
# Users must define their own paths

# datapath<-"user defined path"
# codepath<-"user defined path"

# And users should remove below line that loads in personal paths
source("/Users/kelsey/Github/neon-species-accumulation/configini.R")


library(lubridate)
library(dplyr)
library(vegan)
library(ggplot2)
library(neonUtilities)

#### And users should define their variables for the NEON data download

# data product of interest, below is for Carabids
product<-"DP1.10022.001"

# start and end dates, must be formatted as "YYYY-MM" or NA for all time
start<-"2000-01"
end<-"2023-12"

# Users need to indicate whether they want to load or re-download and format the NEON data 
# These these steps are very time consuming so it is recommended that they are only done if necessary
NewData<-FALSE

# Load in the Data

# Load or download the NEON observation data
if(NewData==TRUE|file.exists(paste0(datapath,"NeonData.Robj"))==FALSE){
  source(paste0(codepath,"DownloadNEONData.R"))
}else{load(file=paste0(datapath,"NeonData.Robj"))}



####
# Pull tables of interest
field<-NeonData3$bet_fielddata
sort<-NeonData3$bet_sorting

# Subset to only traps that were collected
field<-field[which(field$sampleCollected=="Y"),]

# Add year 
field$year<-year(field$collectDate)

# Get the unique site x years list
siteSummary<-data.frame(field %>% group_by(siteID) %>% summarise(years=length(unique(year))))

# Source the randomSpAccum function
source("/Users/kelsey/Github/neon-species-accumulation/code/randomTrapFunction.R")

### iterate over all sites

for (i in 1:nrow(siteSummary)){
  site<-siteSummary$siteID[i]
  field.dat<-field[which(field$siteID==site),]
  sort.dat<-sort[which(sort$siteID==site),]
  sort.dat$year<-year(sort.dat$collectDate)
  trapsPerYear<-field.dat %>% group_by(year) %>% count()
  
  rand.sort.totals<-data.frame(sort.dat %>% group_by(taxonID,year) %>% summarise(individualCount=sum(individualCount)))
  
  # reshape to create typical community matrix
  comm<-reshape(rand.sort.totals[,c("taxonID","year","individualCount")],direction="wide",timevar="taxonID",idvar="year")
  
  # redefine NAs as 0s
  comm[is.na(comm)]<-0
  
  accum<-specaccum(comm)
  
  plot(accum,main=as.character(site))
  fit<-fitspecaccum(accum,"lomolino")
  coef(fit)
  }
  








