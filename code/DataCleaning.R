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

#### Determine incorrect species identifications and propagate correct identifications into main sorting table

# Join expert and parataxonomist tables
para_expert<-left_join(para,expert,join_by("individualID"=="individualID"),suffix=c(".para",".expert"))

# All expert rows should correspond to one and only one para row, but not all para rows will correspond to an expert row
# Investigate duplicates if there is a mismatch between number of rows in para and para_expert
para_expert[which(para_expert$uid.para %in% para_expert$uid.para[which(duplicated(para_expert$uid.para))]),]
# Some duplicates are potentially re-identifications by different experts of the same beetles, leaving all data to be conservative, does not affect overall presence of species at any site

# Subset to instances where expert ID exists to investigate mismatches
para_expert_Exp<-para_expert[which(is.na(para_expert$scientificName.expert)==FALSE),]
mismatch<-para_expert_Exp[which(para_expert_Exp$scientificName.para!=para_expert_Exp$scientificName.expert),]
# see that many mismatches exist, so decide need to propagate expert IDs into main sorting table

# For all pinned specimens choose the ID of either the 
para_expert$chosenID<-para_expert$taxonID.expert
para_expert$chosenID[is.na(para_expert$chosenID)]<-para_expert$taxonID.para[is.na(para_expert$chosenID)]

for (i in 1:nrow()){
  
}

# Put expert ID where exist into sorting table
sort_expert<-left_join(sort,para_expert,join_by("subsampleID"=="subsampleID"),suffix=c(".sort",".para_expert"))



# Examine how often beetles from the same subsample are identified differently be the exper
multiples<-sort_expert$subsampleID[duplicated(sort_expert$subsampleID)]
for (i in 1:length(multiples){
  dupes<-sort_expert[which(sort_expert$subsampleID==multiples[i]),]
  
}


for (i in 1:nrow(expert)){
  # which parataxonomy table beetles are associated with the given pinned beetle identification?
  subID<-para_expert$subsampleID[which(para_expert$individualID==expert$individualID[i])]
  # what is the expert ID for these samples
  sp<-para_expert$taxonID.para[which(para_expert$individualID==expert$individualID[i])]
  # find the related sampleIDs
  sortIDs<-which(sort$subsampleID==subID)# & sort$taxonID==sp)
  print(paste(i,sortIDs))
}


