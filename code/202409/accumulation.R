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
library("iNEXT")

#set seed
set.seed(85705)

# Load in the formatted clean data, or download and create it. 
#Make sure DataCleaning.R is correctly configured

if(NewCleanData==TRUE|file.exists(paste0(datapath,"CleanedData.Robj"))==FALSE){
  source(paste0(codepath,"DataCleaning.R"))}else{load(file=paste0(datapath,"202409_CleanedData.Robj"))}

# Define data of interest
fullData<-CleanedData$fullData
field<-CleanedData$field

# Define the incidence frequency data tables

sites<-unique(field$siteID)

inext<-setNames(vector(mode="list",length=length(sites)),sites)

for (i in 1:length(inext)){
  
  dat<-fullData[which(fullData$siteID==sites[i]),]
  years<-unique(dat$year)
  
  inext[[i]]<-setNames(vector(mode="list",length=length(years)),years)
  inc_freq<-setNames(vector(mode="list",length=length(years)),years)
  
  for (j in 1:length(year)){
    
    datyear<-dat[which(dat$year==years[j]),]
    samps<-unique(datyear$sampleID)
    spp<-unique(datyear$taxonID)
    spp<-spp[is.na(spp)==FALSE]
    inc<-data.frame(matrix(ncol=length(samps),nrow=length(spp)))
    colnames(inc)<-samps
    rownames(inc)<-spp
    
    for (k in 1:nrow(inc)){
      for (l in 1:ncol(inc)){
        inds<-dat$individualCountFinal[which(datyear$sampleID==samps[k] & datyear$taxonID==spp[l])]
        if (length(inds)>0) {inc[k,l]<-sum(inds)}
      }
    }
    
    inc[is.na(inc)]<-0
    
    inext[[i]][[j]]$inc<-inc
    
    presabs<-inc
    presabs[presabs>1]<-1
    
    input<-c(ncol(presabs),as.vector(rowSums(presabs)))
    inext[[i]][[j]]$inc_freq<-input
    inc_freq[[j]]<-input
  }
    

    out<-iNEXT(inc_freq,datatype='incidence_freq',knot=20,endpoint=ncol(presabs)*3)
    inext[[i]]$out<-out
    
  }
  
}





