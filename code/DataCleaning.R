### This is the primary file for creating tables of cleaned data for analysis

#### BEFORE RUNNING
# Users must define their own paths

# datapath<-"user defined path"
# codepath<-"user defined path"

#And users must have downloaded the data file as desired for their work

# And users should remove below line that loads in personal paths
source("/Users/kelsey/Github/neon-species-accumulation/configini.R")

# load NEON carabid taxonomy table
taxa<-read.csv(paste0(datapath,"CarabidTaxonomicList_March2024.csv"))

# Load or download the NEON observation data
source(paste0(codepath,"DownloadNEONData.R"))

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
for (i in 1:nrow(pinned)){
  row<-which(main$subsampleID==pinned$subsampleID[i])
  main$individualCount[row]<-main$individualCount[row]-1
  newrecord<-main[row,]
  newrecord$individualCount<-1
  newrecord$taxonID<-pinned$taxonID[i]
  main<-rbind(main,newrecord)
}

# Drop all instances in which identification was not Carabid or was only to the family level

