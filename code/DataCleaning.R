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

# Create records for the pinned individuals in the main table
#### This is messing it up because it's renaming all of the records
for (i in 1:nrow(pinned)){
  row<-which(main$subsampleID==pinned$subsampleID[i])
  main$individualCount[row]<-main$individualCount[row]-1
  newrecords<-main[row,]
  newrecords$individualCount<-1
  newrecords$taxonID<-pinned$taxonID[i]
  main<-rbind(main,newrecords)
  print(paste(i,nrow(newrecords)))
}

# Replace subspecies with species level ID in the taxa table
taxa$sciName<-taxa$scientificName
taxa$sciName[which(taxa$taxonRank=="subspecies")]<-word(taxa$sciName, 1,2, sep=" ")

# Join main and taxonomy tables & Subset taxonomy table to fields of interest, removes subspecies
fullData<-left_join(main,
                taxa[,c("taxonID","family","sciName")],
                join_by("taxonID"=="taxonID"))

# Drop all instances in which identification was not Carabid or was only to the family level
fullData<-fullData[which(fullData$family=="Carabidae"),]
fullData<-fullData[-which(fullData$sciName %in% ("Carabidae sp.","Carabidae spp.")),]

# Sampling effort

sort<-sort[which(sort$sampleType %in% c("carabid","other carabid")),]
sort$year<-format(as.Date(sort$collectDate),'%Y')
sort<-left_join(sort,field,join_by("sampleID"=="sampleID"))
indivsEffort<-data.frame(sort %>% 
                           group_by(siteID.x,nlcdClass,year.x) %>% 
                           summarise(indivs=sum(individualCount)))



#Make a table of trapping days x siteHabitat x year
fullData$year<-format(as.Date(fullData$collectDate),'%Y')
dhy_effort<-fullData[which(duplicated(fullData$sampleID)==FALSE),]

