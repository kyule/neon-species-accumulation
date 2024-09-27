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
library("codyn")
library("reshape")
library("rarestR")


#set seed
set.seed(85705)

# Load in the formatted clean data, or download and create it. 
#Make sure DataCleaning.R is correctly configured

if(NewCleanData==TRUE|file.exists(paste0(datapath,"FullAndCleanData.Robj"))==FALSE){
  source(paste0(codepath,"DataCleaningV2.R"))}else{load(file=paste0(datapath,"FullAndCleanData.Robj"))}

# Pull data of interest
fullData<-FullAndCleanData$fullData
fullData$year<-year(fullData$collectDate)

# Function for creating incidence data frame

inc_data<-function(input.data){
  
  # find the unique samples and spp for that year
  samps<-unique(input.data$sampleID)
  spp<-unique(input.data$sciName)
  spp<-spp[is.na(spp)==FALSE]
  
  return.df<-list()
  
  return.df$traps<-samps
  return.df$spp<-spp
  
  # initiate an incidence data frame
  inc<-data.frame(matrix(ncol=length(samps),nrow=length(spp)))
  colnames(inc)<-samps
  rownames(inc)<-spp
  inc[is.na(inc)]<-0
  
  # Input individuals into incidence data frame
  for (k in 1:nrow(input.data)){
    row<-which(rownames(inc)==input.data$sciName[k])
    if (length(row)>0){
      col<-which(colnames(inc)==input.data$sampleID[k])
      if (length(col)>0){
        inc[row,col]<-inc[row,col]+input.data$finalIndivCount[k]
      }
    }
  }
  
  
  return.df$inc<-inc
  
  # Convert to presence absence
  presabs<-inc
  presabs[presabs>1]<-1
  
  # Create and save the incidence frequency input data for iNext
  input<-c(ncol(presabs),as.vector(rowSums(presabs)))
  return.df$inc_freq<-input
  
  return(return.df)
  
}


# Create overarching data frame for analysis

sites<-unique(fullData$siteID)
results<-setNames(vector(mode="list",length=length(sites)),sites)


# loop through the field sites
for (i in 1:length(results)){
  
  print(paste0("START: ", sites[i]))
  
  #subset data to the field site of interest and determine the years of analysis
  dat<-fullData[which(fullData$siteID==sites[i]),]
  years<-unique(dat$year)
  
  # Create list structures for the site
  results[[i]]<-setNames(vector(mode="list",length=(length(years)+1)),c(years,"full"))
  inc_freq<-setNames(vector(mode="list",length=(length(years)+1)),c(years,"full"))
  
  
  for (j in 1:length(years)){
    # loop through the years of sampling at the site
    print(paste(sites[i],years[j]),sep=": ")
    
    # pull out the data for the specific year
    datyear<-dat[which(dat$year==years[j]),]
    
    # results by year
    results[[i]][[j]]<-inc_data(datyear)
    inc_freq[[j]]<-results[[i]][[j]]$inc_freq
    
  }
  
  # Now do the same across the full data for the site
    print(paste0(sites[i]," Full Data"))
    
    results[[i]][[j+1]]<-inc_data(dat)
    inc_freq[[j+1]]<-results[[i]][[j+1]]$inc_freq
    
    # Conduct analysis on the communities, each year + the full data are treated as different 'assemblages'
    print(paste0("****iNEXT analysis: ",sites[i]))
    out<-try({
      iNEXT(inc_freq,q=c(0,1,2),datatype='incidence_freq',knot=20,endpoint=ncol(presabs)*3)},
      silent=TRUE)
    if (inherits(out, "try-error")) {
      results[[i]]$out <- paste("Error")
    } else {
      results[[i]]$out<-out
    }
    
    
    ### Calculate turnover measures
    
    print(paste0("turnover: ",sites[i]))
    turnover.df<-data.frame(dat %>% group_by(year,sciName) %>% summarise(count=sum(finalIndivCount)))
    turnover.df<-turnover.df[turnover.df$count>=1,]
    total<-turnover(turnover.df,time.var="year",species.var="sciName",abundance.var="count",metric="total")
    appear<-turnover(turnover.df,time.var="year",species.var="sciName",abundance.var="count",metric="appearance")
    disappear<-turnover(turnover.df,time.var="year",species.var="sciName",abundance.var="count",metric="disappearance")
    turnover<-left_join(total,appear,join_by("year"=="year"))
    turnover<-left_join(turnover,disappear,join_by("year"=="year"))
    
    results[[i]]$turnover<-turnover
    
    # dissimilarity measures
    
    print(paste0("dissimilarity: ",sites[i]))
    dat$siteyear<-paste(dat$siteID,dat$year,sep="_")
    dis.dat<-data.frame(dat %>% group_by(siteyear,sciName) %>% summarise(n=sum(finalIndivCount)))
    dis.dat<-dis.dat[which(dis.dat$n>0),]
    dis.dat<-acast(dis.dat, siteyear~sciName, value.var="n")
    dissim<-ess(dis.dat)
    
    results[[i]]$dissim<-dissim
    
}



  





