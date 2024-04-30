### This is the primary file for creating tables of cleaned data for analysis

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

# load NEON carabid taxonomy table
taxa<-read.csv(paste0(datapath,"CarabidTaxonomicList_March2024.csv"))

# Load or download the NEON observation data
if(NewData==TRUE|file.exists(paste0(datapath,"NeonData.Robj"))==FALSE){
  source(paste0(codepath,"DownloadNEONData.R"))
  }else{load(file=paste0(datapath,"NeonData.Robj"))}

# Pull tables of interest
field<-NeonData$bet_fielddata
sort<-NeonData$bet_sorting
para<-NeonData$bet_parataxonomistID
expert<-NeonData$bet_expertTaxonomistIDProcessed

#  Subset field to whether sample was collected to start main analysis file
main<-field[which(field$sampleCollected=="Y"),]

#subset to columnns of interest
main<-main[,c("plotID","siteID","nlcdClass","setDate","collectDate","sampleID")]

#Subset sort table to carabids, columns of interest and join to field table
main<-left_join(main,
                sort[which(sort$sampleType %in% c("other carabid","carabid")),c("sampleID","subsampleID","taxonID","individualCount")],
                join_by("sampleID"=="sampleID"))

#### Determine incorrect species identifications and propagate correct identifications into main sorting table

# Join necessary columns of expert and parataxonomist tables
pinned<-left_join(para[,c("subsampleID","individualID","taxonID")],
                  expert[c("individualID","taxonID")],
                  join_by("individualID"=="individualID"),
                  suffix=c(".para",".expert"))

# For all pinned specimens choose the ID of the expert when exists, parataxonomist otherwise
pinned$taxonID<-pinned$taxonID.expert
pinned$taxonID[is.na(pinned$taxonID)]<-pinned$taxonID.para[is.na(pinned$taxonID)]

### Create records for the pinned individuals in the main table
# choose which subsamples have pinned specimens associated with them
subsToUpdate<-unique(pinned$subsampleID)
main$individualCountFinal<-main$individualCount
  
for (i in 1:length(subsToUpdate)){
  row<-which(main$subsampleID==subsToUpdate[i])
  pinnedIndivs<-pinned[which(pinned$subsampleID==subsToUpdate[i]),]
  main$individualCountFinal[row]<-main$individualCountFinal[row]-nrow(pinnedIndivs)
  newrecords<-data.frame(plotID=main$plotID[row],siteID=main$siteID[row],
                         nlcdClass=main$nlcdClass[row],setDate=main$setDate[row],
                         collectDate=main$collectDate[row],sampleID=main$sampleID[row],
                         subsampleID=main$subsampleID[row],taxonID=pinnedIndivs$taxonID,
                         individualCount=1,individualCountFinal=1)
  main<-rbind(main,newrecords)
}

# Replace NAs in individualCount with 0s
fullData<-main[ , !(names(main) %in% "individualCount")] # replace to not risk losing what took a long time to run :)
fullData$individualCountFinal[is.na(fullData$individualCountFinal)]<-0

# Replace subspecies with species level ID in the taxa table
taxa$sciName<-taxa$scientificName
taxa$sciName[which(taxa$taxonRank=="subspecies")]<-word(taxa$sciName[which(taxa$taxonRank=="subspecies")], 1,2, sep=" ")
taxa$sciName<-str_replace(taxa$sciName," sp.","")
taxa$sciName<-str_replace(taxa$sciName," spp.","")

# Join main and taxonomy tables & Subset taxonomy table to fields of interest, removes subspecies
fullData<-left_join(fullData,
                taxa[,c("taxonID","family","sciName")],
                join_by("taxonID"=="taxonID"))

# Drop all instances in which identification was not Carabid or was only to the family level
fullData<-fullData[which(fullData$family %in% c(NA,"Carabidae")),]
fullData<-fullData[-which(fullData$sciName %in% c("Carabidae","Carabidae")),]


##### Sampling effort by year,site,nlcdclass

# Table by number of individuals
fullData$year<-format(as.Date(fullData$collectDate),'%Y')
effort<-data.frame(fullData %>% 
                           group_by(siteID,nlcdClass,year) %>% 
                           summarise(indivs=sum(individualCountFinal)))

# Table of trapping days x siteHabitat x year
dhy_effort<-fullData[which(duplicated(fullData$sampleID)==FALSE),]
dhy_effort$trapDay<-as.numeric(as.Date(dhy_effort$collectDate) - as.Date(dhy_effort$setDate))
daysEffort<-data.frame(dhy_effort %>%
                         group_by(siteID,nlcdClass,year) %>%
                         summarise(days=sum(trapDay)))

# Get individuals and days together
effort<-full_join(effort,daysEffort,join_by("siteID"=="siteID","nlcdClass"=="nlcdClass","year"=="year"))

### Save these data for analyses
CleanedData<-list(fullData=fullData,effort=effort)
save(CleanedData,file=paste0(datapath,"CleanedData.Robj"))

#Remove unnecessary data

rm(list=ls()[! ls() %in% c("CleanedData")])
