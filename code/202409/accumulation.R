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
library("stringr")
library("iNEXT")


#set seed
set.seed(85705)

# Load in the formatted clean data, or download and create it. 
#Make sure DataCleaning.R is correctly configured

if(NewCleanData==TRUE|file.exists(paste0(datapath,"CleanedData.Robj"))==FALSE){
  source(paste0(codepath,"DataCleaningV2.R"))}else{load(file=paste0(datapath,"FullAndCleanData.Robj"))}

# Pull data of interest
fullData<-FullAndCleanData$fullData
fullData$year<-year(fullData$collectDate)

# Create overarching data frame for analysis

sites<-unique(fullData$siteID)
inext<-setNames(vector(mode="list",length=length(sites)),sites)


# loop through the field sites
for (i in 1:length(inext)){
  
  print(paste0("START: ", sites[i]))
  
  #subset data to the field site of interest and determine the years of analysis
  dat<-fullData[which(fullData$siteID==sites[i]),]
  years<-unique(dat$year)
  
  # Create list structures for the site
  inext[[i]]<-setNames(vector(mode="list",length=(length(years)+1)),c(years,"full"))
  inc_freq<-setNames(vector(mode="list",length=(length(years)+1)),c(years,"full"))
  
  # loop through the years of sampling at the site
  for (j in 1:length(years)){
    print(paste(sites[i],years[j]),sep=": ")
    
    # pull out the data for the specific year
    datyear<-dat[which(dat$year==years[j]),]
    
    # find the unique samples and spp for that year
    samps<-unique(datyear$sampleID)
    spp<-unique(datyear$sciName)
    spp<-spp[is.na(spp)==FALSE]
    
    inext[[i]][[j]]$traps<-samps
    inext[[i]][[j]]$spp<-spp
    
    # initiate an incidence data frame
    inc<-data.frame(matrix(ncol=length(samps),nrow=length(spp)))
    colnames(inc)<-samps
    rownames(inc)<-spp
    
    # Input individuals into incidence data frame
    for (k in 1:nrow(inc)){
      for (l in 1:ncol(inc)){
        inds<-dat$finalIndivCount[which(datyear$sampleID==samps[k] & datyear$sciName==spp[l])]
        if (length(inds)>0) {inc[k,l]<-sum(inds)}
      }
    }
    
    inc[is.na(inc)]<-0
    inext[[i]][[j]]$inc<-inc
    
    # Convert to presence absence
    presabs<-inc
    presabs[presabs>1]<-1
    
    # Create and save the incidence frequency input data for iNext
    input<-c(ncol(presabs),as.vector(rowSums(presabs)))
    inext[[i]][[j]]$inc_freq<-input
    inc_freq[[j]]<-input
  }
  
  print(paste0(sites[i]," Full Data"))
    samps<-unique(dat$sampleID)
    spp<-unique(dat$taxonID)
    spp<-spp[is.na(spp)==FALSE]
    inc<-data.frame(matrix(ncol=length(samps),nrow=length(spp)))
    colnames(inc)<-samps
    rownames(inc)<-spp
    
    for (k in 1:nrow(inc)){
      for (l in 1:ncol(inc)){
        inds<-dat$individualCountFinal[which(dat$sampleID==samps[k] & dat$taxonID==spp[l])]
        if (length(inds)>0) {inc[k,l]<-sum(inds)}
      }
    }
    
    inc[is.na(inc)]<-0
    
    inext[[i]][[j+1]]$inc<-inc
    
    presabs<-inc
    presabs[presabs>1]<-1
    
    input<-c(ncol(presabs),as.vector(rowSums(presabs)))
    inext[[i]][[j+1]]$inc_freq<-input
    inc_freq[[j+1]]<-input
    print(paste0("iNEXT ",sites[i]))
    out<-iNEXT(inc_freq,q=c(0,1,2),datatype='incidence_freq',knot=20,endpoint=ncol(presabs)*3)
    inext[[i]]$out<-out
    
}


  





