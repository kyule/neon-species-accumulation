### This is the primary file for analyzing NEON species accumulation data

#### BEFORE RUNNING
# Users must define their own paths

# datapath<-"user defined path"
# codepath<-"user defined path"

# Users need to indicate whether they want to load or re-download and format the NEON data 
# These these steps are very time consuming so it is recommended that they are only done if necessary
NewCleanData<-FALSE

# And users should remove below line that loads in personal paths
source("/Users/kelsey/Github/neon-species-accumulation/configini.R")

# Load some necessary packages
library("dplyr")
library("vegan")

# Load in the formatted clean data, or download and create it. 
#Make sure DataCleaning.R is correctly configured

if(NewCleanData==TRUE|file.exists(paste0(datapath,"CleanedData.Robj"))==FALSE){
  source(paste0(codepath,"DataCleaning.R"))}else{load(file=paste0(datapath,"CleanedData.Robj"))}

# Define data of interest
fullData<-CleanedData$fullData
effort<-CleanedData$effort

# Creat function to reformat data for community matrix input into vegan package ### SOMETHING IS WRONG
createComm<-function(data,site,nlcd){
  community<-data[which(data$siteID==site & data$nlcdClass==nlcd),]
  community<-community[!is.na(community$taxonID),]
  community$comm<-paste(community$siteID,community$nlcdClass,community$year,sep=".")
  community<-community[,colnames(community) %in% c("comm","sciName","individualCountFinal")]
  communityMatrix<-reshape(community,direction="wide",timevar = "sciName",idvar="comm")
  colnames(communityMatrix)<-str_replace(colnames(communityMatrix),"individualCountFinal.","")
  communityMatrix[is.na(communityMatrix)]<-0
  return(communityMatrix)
}
