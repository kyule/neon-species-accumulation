
### This is the primary file for creating tables of cleaned data for analysis

#### BEFORE RUNNING
# Users must define their own paths

# datapath<-"user defined path"
# codepath<-"user defined path"

#And users must have downloaded the data file as desired for their work

# And users should remove below line that loads in personal paths
source("/Users/kelsey/Github/neon-species-accumulation/configini.R")

# Load or download the NEON observation data and clean it
source(paste0(codepath,"DownloadNEONData.R"))
source(paste0(codepath,"DataCleaning.R"))

# Load necessary package
library("dplyr")

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

