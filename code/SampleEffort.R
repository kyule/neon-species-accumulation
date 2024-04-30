### This is the primary file for determining sampling effort

#### BEFORE RUNNING
# Users must define their own paths

# datapath<-"user defined path"
# codepath<-"user defined path"

#And users must have configured the DownloadNEONData.R file as desired for their work
# Indicate whether new data should be downloaded, downloads are time consuming so this is not recommended
NewData<-FALSE

# And users should remove below line that loads in personal paths
source("/Users/kelsey/Github/neon-species-accumulation/configini.R")

# load necessary libraries
library(dplyr)
library(stringr)

# Load or download the NEON observation data
if(NewData==TRUE|file.exists(paste0(datapath,"NeonData.Robj"))==FALSE){
  source(paste0(codepath,"DownloadNEONData.R"))
}else{load(file=paste0(datapath,"NeonData.Robj"))}

# Field sampling effort in trap nights, site x habitat x year
field<-field[which(field$sampleCollected=="Y"),]
field$year<-format(as.Date(field$collectDate),'%Y')
daysEffort<-data.frame(field %>% 
                         group_by(siteID,nlcdClass,year) %>% 
                         summarise(days=sum(trappingDays)))

# Pull tables of interest
field<-NeonData$bet_fielddata
sort<-NeonData$bet_sorting

# Sampling effort in individuals, site x habitat x year
sort<-sort[which(sort$sampleType %in% c("carabid","other carabid")),]
sort$year<-format(as.Date(sort$collectDate),'%Y')
sort<-left_join(sort,field,join_by("sampleID"=="sampleID"))
indivsEffort<-data.frame(sort %>% 
                           group_by(siteID.x,nlcdClass,year.x) %>% 
                           summarise(indivs=sum(individualCount)))

