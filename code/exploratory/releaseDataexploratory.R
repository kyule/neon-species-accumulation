
### This is the primary file for creating site based results necessary for NEON species accumulation work

#### BEFORE RUNNING
# Users must define their own paths

# datapath<-"user defined path"
# codepath<-"user defined path"

# And users should remove below line that loads in personal paths
source("/Users/kelsey/Github/neon-species-accumulation/configini.R")

# Load some necessary packages
library("dplyr")
library("stringr")
library("iNEXT")
library("codyn")
library("reshape2")

#set seed
set.seed(85705)

# set completeness threshold
comp.thresh<-0.1

# Load in the formatted clean data

load(file=paste0(datapath,"CleanNEONData.Robj"))

# Pull data of interest
fullData<-FullAndCleanData$fullData
completeness<-FullAndCleanData$completeness

# Function for creating incidence data frame

inc_data<-function(input.data){
  
  # find the unique samples and spp
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
  
  #iNext results
  
  warnings_list <- character()  # Store warnings
  
  out <- try({
    withCallingHandlers(
      iNEXT(input, q = c(0, 1), datatype = 'incidence_freq', knot = 30, endpoint = ncol(presabs) * 10),
      warning = function(w) {
        warnings_list <<- c(warnings_list, conditionMessage(w))  # Capture warning messages
        invokeRestart("muffleWarning")  # Prevent warning from displaying
      }
    )
  }, silent = TRUE)
  
  # Store results, capturing both errors and warnings
  if (inherits(out, "try-error")) {
    return.df$out <- "Error"
    return.df$warning <- NA  # No warning if there was an error
  } else {
    return.df$out <- out
    return.df$warning <- if (length(warnings_list) > 0) paste(warnings_list, collapse = "; ") else NA
  }
  
  
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
  rems<-completeness$year[which(completeness$siteID==sites[i] & completeness$propRem>comp.thresh)]
  rems<-c(rems,completeness$year[which(completeness$siteID==sites[i] & is.nan(completeness$propRem))])
  # choose to remove years for which more than the completeness threshhold of the beetles were not carabids and/or not identified to a lower resolution than "carabid"
  if (length(rems)>0) {dat<-dat[-which(dat$year %in% rems),]}
  years<-unique(dat$year)
  
  # Create list structures for the site
  results[[i]]<-setNames(vector(mode="list",length=(length(years)+1)),c(years,"full"))
  inc_freq<-setNames(vector(mode="list",length=(length(years)+1)),c(years,"full"))
  
  
  for (j in 1:length(years)){
    # loop through the years of sampling at the site
    print(paste(sites[i],years[j]),sep=": ")
    
    # pull out the data for the specific year
    datyear<-dat[which(dat$year<=years[j]),]
    
    # results by year
    results[[i]][[j]]<-inc_data(datyear)
    inc_freq[[j]]<-results[[i]][[j]]$inc_freq
    
  }
  
  # Now do the same across the full data for the site
  print(paste0(sites[i]," Full Data"))
  
  results[[i]][[j+1]]<-inc_data(dat)
  inc_freq[[j+1]]<-results[[i]][[j+1]]$inc_freq
  
  
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
  
  
}

save(results,file=paste0(datapath,"iNEXTandTurnoverResults.Robj"))



##### explore warnings

warnings_list <- lapply(results, function(site) {
  lapply(site, function(year) {
    if (!is.null(year$warning)) return(year$warning) else return(NULL)
  })
})

# Flatten the list to remove NULL values
warnings_flat <- data.frame(unlist(warnings_list, recursive = TRUE, use.names = TRUE))

# View warnings
warnings_flat<-data.frame(site=row.names(warnings_flat),warning=warnings_flat)
warnings_flat$year<-sapply(strsplit(row.names(warnings_flat),"[.]"),"[",2)
warnings_flat$site<-sapply(strsplit(row.names(warnings_flat),"[.]"),"[",1)
