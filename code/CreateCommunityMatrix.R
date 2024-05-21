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
library("stringr")

#set seed
set.seed(85705)

# Load in the formatted clean data, or download and create it. 
#Make sure DataCleaning.R is correctly configured

if(NewCleanData==TRUE|file.exists(paste0(datapath,"CleanedData.Robj"))==FALSE){
  source(paste0(codepath,"DataCleaning.R"))}else{load(file=paste0(datapath,"CleanedData.Robj"))}

# Define data of interest
fullData<-CleanedData$fullData
effort<-CleanedData$effort

# Create function to reformat data for community matrix input into vegan package ### SOMETHING IS WRONG
createComm<-function(data,effort,site,nlcd,standardized){
  community<-data[which(data$siteID==site & data$nlcdClass==nlcd),]
  community<-community[community$individualCountFinal!=0,]
  community$comm<-paste(community$siteID,community$nlcdClass,community$year,sep=".")
  community<-community[,colnames(community) %in% c("comm","sciName","individualCountFinal")]
  community.summary<-data.frame(community %>% group_by(comm,sciName) %>% summarise(individualCountFinal=sum(individualCountFinal)))
  communityMatrix<-reshape(community.summary,direction="wide",timevar = "sciName",idvar="comm")
  colnames(communityMatrix)<-str_replace(colnames(communityMatrix),"individualCountFinal.","")
  communityMatrix[is.na(communityMatrix)]<-0
  if(standardized=="years"){
    row.names(communityMatrix)<-sapply(strsplit(communityMatrix$comm,"[.]"),'[',3)
    return(communityMatrix[,-1])}
  if(standardized=="days"){
    eff<-effort[which(effort$siteID==site & effort$nlcdClass==nlcd),]
    eff$comm<-paste(eff$siteID,eff$nlcdClass,eff$year,sep=".")
    communityMatrix<-left_join(communityMatrix,eff[,c("comm","days")],join_by("comm"=="comm"))
    communityMatrix[,-c(1,ncol(communityMatrix))]<-communityMatrix[,-c(1,ncol(communityMatrix))]/communityMatrix$days
    row.names(communityMatrix)<-sapply(strsplit(communityMatrix$comm,"[.]"),'[',3)
    return(communityMatrix[,-c(1,ncol(communityMatrix))])
  }
  if(standardized=="individuals"){
    eff<-effort[which(effort$siteID==site & effort$nlcdClass==nlcd),]
    eff$comm<-paste(eff$siteID,eff$nlcdClass,eff$year,sep=".")
    communityMatrix<-left_join(communityMatrix,eff[,c("comm","indivs")],join_by("comm"=="comm"))
    communityMatrix[,-c(1,ncol(communityMatrix))]<-communityMatrix[,-c(1,ncol(communityMatrix))]/communityMatrix$indivs
    row.names(communityMatrix)<-sapply(strsplit(communityMatrix$comm,"[.]"),'[',3)
    return(communityMatrix[,-c(1,ncol(communityMatrix))])
  }
}

a<-createComm(fullData,effort,"SRER","shrubScrub","days")
b<-specaccum(a)
plot(b)

c<-createComm(fullData,effort,"SRER","shrubScrub","years")
d<-specaccum(c)
plot(d)

e<-createComm(fullData,effort,"SRER","shrubScrub","individuals")
f<-specaccum(e)
plot(f)

raremax<-min(rowSums(c))
rarefy(c,raremax)
rarecurve(c,1,100)

metaMDS(c)->ord
plot(ord)


#### Just for fun, all communities in an ordination
fullComms<-data.frame(fullData %>% group_by(siteID,sciName) %>% summarise(n=sum(individualCountFinal)))
### Something is wrong because individual counts of 0 shouldn't be showing up
communityMatrix<-reshape(fullComms,direction="wide",timevar = "sciName",idvar="siteID")
communityMatrix[is.na(communityMatrix)]<-0
row.names(communityMatrix)<-communityMatrix[,1]
metaMDS(communityMatrix[,-1])->ord
plot(ord,display='sites',type='t')
