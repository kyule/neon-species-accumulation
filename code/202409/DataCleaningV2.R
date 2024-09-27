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

#Subset sort table to carabids, columns of interest and join to field table
main<-full_join(main,
                sort[which(sort$sampleType %in% c("other carabid","carabid")),c("sampleID","subsampleID","taxonID","individualCount")],
                join_by("sampleID"=="sampleID"))

# Replace subspecies with species level ID in the taxa table
taxa$sciName<-taxa$scientificName
taxa$sciName[which(taxa$taxonRank=="subspecies")]<-word(taxa$sciName[which(taxa$taxonRank=="subspecies")], 1,2, sep=" ")
taxa$sciName<-str_replace(taxa$sciName," sp.","")
taxa$sciName<-str_replace(taxa$sciName," spp.","")
taxa$sciName[which(taxa$taxonRank=="subgenus")]<-word(taxa$sciName[which(taxa$taxonRank=="subgenus")], 1,1, sep=" ")

# Replace taxonID with sciName used in analysis in all tables
main<-left_join(main,taxa[,c("taxonID","sciName")],join_by("taxonID"=="taxonID"))
main<-subset(main,select=-taxonID)
sort<-left_join(sort,taxa[,c("taxonID","sciName")],join_by("taxonID"=="taxonID"))
sort<-subset(sort,select=-taxonID)
para<-left_join(para,taxa[,c("taxonID","sciName")],join_by("taxonID"=="taxonID"))
para<-subset(para,select=-taxonID)
expert<-left_join(expert,taxa[,c("taxonID","sciName")],join_by("taxonID"=="taxonID"))
expert<-subset(expert,select=-taxonID)

#### Determine incorrect species identifications and propagate correct identifications into main sorting table

# Join necessary columns of expert and parataxonomist tables
pinned<-left_join(expert[c("individualID","sciName")],
                  para[,c("subsampleID","individualID","sciName")],
                  join_by("individualID"=="individualID"),
                  suffix=c(".expert",".para"))

#### Determine misidentifications

misIDs<-pinned[which(pinned$sciName.expert!=pinned$sciName.para),]
nrow(misIDs)/nrow(pinned) # total misID rate for all expert ID beetles is ~19%

#### Propagate correct identifications into main sorting table
# Note that correcting only the individuals in the sorting table that are explicitly known to be wrong is most conservative in that it maximizes the number of species found within single trap

# choose which subsamples have pinned specimens associated with misidentified individuals and define indiv counts
misIDsubs<-data.frame(misIDs %>% group_by(subsampleID,sciName.expert,sciName.para) %>% count())
main$individualCount[is.na(main$individualCount)]<-0
main$finalIndivCount<-main$individualCount



# Add new records for misidentified subsample individuals

fullData<-main
for (i in 1:nrow(misIDsubs)){
  print(i)
  # find the matching row
  row<-which(fullData$subsampleID==misIDsubs$subsampleID[i] & fullData$sciName==misIDsubs$sciName.para[i])
  # update the main records individual count
  if(length(row)==1){
    fullData$finalIndivCount[row]<-fullData$finalIndivCount[row]-misIDsubs$n[i]
  # create new records
    newrecords<-fullData[row,]
    newrecords$sciName<-misIDsubs$sciName.expert[i]
    newrecords$finalIndivCount<-misIDsubs$n[i]
  #add new records
    fullData<-rbind(fullData,newrecords)
  }
}



# Identify sciNames to remove

removeTaxa<-taxa$sciName[which(taxa$family!="Carabidae")]
removeTaxa<-c(removeTaxa,taxa$sciName[grep("Carabidae",taxa$sciName)])

# Determine number of traps by year, proportion of individuals not sufficiently ID'd, and error rate by site x year combo

field<-field[which(field$sampleCollected=="Y"),]
field$year<-year(field$collectDate)
completeness<-data.frame(field %>% group_by(siteID,year) %>% summarise(traps=length(unique(sampleID))))

fullData$year<-year(fullData$collectDate)
sum.fullTaxa<-data.frame(fullData %>% group_by(siteID,year) %>% summarise(all=sum(finalIndivCount),rems=sum(finalIndivCount[which(sciName %in% removeTaxa)])))
sum.fullTaxa$propRem<-sum.fullTaxa$rems/sum.fullTaxa$all

completeness<-full_join(completeness,sum.fullTaxa,join_by("siteID"=="siteID","year"=="year"))

# Replace sciName for records identified only to family or which are not carabids with NA
fullData$sciName[which(fullData$sciName %in% removeTaxa)]<-NA

#save all of the data for other use
FullAndCleanData<-list(fullData=fullData,field=field,sort=sort,expert=expert,para=para,taxa=taxa)
save(FullAndCleanData,file=paste0(datapath,"FullAndCleanData.Robj"))







