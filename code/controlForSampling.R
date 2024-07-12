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
end<-"2024-05"

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
field<-NeonData$bet_fielddata
sort<-NeonData$bet_sorting

# Subset to only traps that were collected
field<-field[which(field$sampleCollected=="Y"),]

# Add year 
field$year<-year(field$collectDate)

# Get the unique site x years list
siteSummary<-data.frame(field %>% group_by(siteID) %>% summarise(years=length(unique(year))))
siteSummary$richness<-NA
siteSummary$obs.asym<-NA
siteSummary$obs.mid<-NA
siteSummary$mean25.asym<-NA
siteSummary$var25.asym<-NA
siteSummary$mean25.mid<-NA
siteSummary$var25.mid<-NA

# Source the randomSpAccum function
source("/Users/kelsey/Github/neon-species-accumulation/code/randomTrapFunction.R")

### iterate over all sites

for (i in 1:nrow(siteSummary)){
  site<-siteSummary$siteID[i]
  field.dat<-field[which(field$siteID==site),]
  sort.dat<-sort[which(sort$siteID==site),]
  trapsPerYear<-field.dat %>% group_by(year) %>% count()
  #n.25<-randomSpAccum(field.dat,trapsPerYear,sort.dat,25,1)
  #n.100<-randomSpAccum(field.dat,trapsPerYear,sort.dat,100,100)
  #max<-randomSpAccum(field.dat,trapsPerYear,sort.dat,min(trapsPerYear$n),100)
  obs<-randomSpAccum(field.dat,trapsPerYear,sort.dat,max(trapsPerYear$n),1)
  
  siteSummary$richness[i]<-length(unique(sort.dat$taxonID[which(sort.dat$sampleType=="carabid")]))
  siteSummary$obs.asym[i]<-obs$asym[1]
  siteSummary$obs.mid[i]<-obs$mid[1]
  
  ggplot(data=obs[1:length(unique(obs$year)),],aes(x=year,y=richness,group=iter,color='red'))+
    geom_line(linewidth=1.5)+
    title(site)+
    theme_minimal()+
    theme(legend.position = 'none')+ 
    ylim(0, (max(obs$richness)+10))
    #+ geom_line(data=n.25,aes(x=year,y=richness,group=iter),color='lightgray',linewidth=0.25)+
    #geom_line(data=n.100,aes(x=year,y=richness,group=iter),color='darkgray',linewidth=0.25)+
    #geom_line(data=max,aes(x=year,y=richness,group=iter),color='blue',linewidth=0.25)
  
}








