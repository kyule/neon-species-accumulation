# Function for running iterations using randomly selected traps

randomSpAccum<-function(field.dat,traps.peryear,sort.dat,n.traps,n.iter){
  
  set.seed(85705)
  
  # subset to carabids
  sort.dat<-sort.dat[which(sort.dat$sampleType=="carabid"),]
  
  # create data frame for species accumulation work
  accum.df<-data.frame(iter=c(),n.traps=c(),year=c(),richness=c(),asym=c(),mid=c())
  
  # cycle through the number of iterations
  for (i in 1:n.iter){
    
    #randomly reorder field data
    field.iter<-field.dat[sample(1:nrow(field.dat),nrow(field.dat)),]
    
    # create data frame for selection of traps for this iteration
    rand.field<-data.frame()
    
    # cycle through the number of years
    for (j in 1:nrow(traps.peryear)){
      #define the year
      year<-traps.peryear$year[j]
      #select traps for the year
      traps<-field.iter[which(field.iter$year==year),]
      #subset to number of traps being sampled which is n.traps unless n.traps is larger than the total number of traps
      rand.traps<-traps[1:min(n.traps,nrow(traps)),]
      # add to the data the randomized data frame
      rand.field<-rbind(rand.field,rand.traps)
    }
    
    # find the years of sampling
    field.years<-unique(rand.field$year)
    # select the sorting table values based on the field samples chosen
    rand.sort<-sort.dat[which(sort.dat$sampleID %in% rand.field$sampleID),]
    #add year to the sorting table
    rand.sort$year<-year(rand.sort$collectDate)
    # are any years missing? (i.e. no carabids were caught)
    miss<-field.years[-which(field.years %in% unique(rand.sort$year))]
    # add empty data for the missing years
    if(length(miss)>0){
      #create skeletal data records for the missing years
      empty.data<-rand.sort[1:length(miss),]
      # set individual count to 0 to indicate that nothing was caught
      empty.data$individualCount<-0
      # update the years to the missing years
      empty.data$year<-miss
      # add it to the rand.sort table
      rand.sort<-rbind(rand.sort,empty.data)
    }
    
    # Count the number of traps by year
    count.traps<-data.frame(rand.field %>% group_by(year) %>% count())$n
    
    # remove the NA specimen counts
    rand.sort<-rand.sort[!is.na(rand.sort$individualCount),]
    
    # summarize by taxon and year and get the total abundance
    rand.sort.totals<-data.frame(rand.sort %>% group_by(taxonID,year) %>% summarise(individualCount=sum(individualCount)))
    
    # reshape to create typical community matrix
    comm<-reshape(rand.sort.totals[,c("taxonID","year","individualCount")],direction="wide",timevar="taxonID",idvar="year")
    
    # redefine NAs as 0s
    comm[is.na(comm)]<-0

    #  species accumulation from vegan package
    accum<-specaccum(comm)
    fit<-fitspecaccum(accum,"lomolino")
    asym<-coef(fit)[1]
    mid<-coef(fit)[2]

    # create a data frame of results for this iteration
    df<-data.frame(iter=rep(i,nrow(traps.peryear)),n.traps=count.traps,year=accum$sites,richness=accum$richness,asym=asym,mid=mid)
    
    # add this iterations results to the overall results
    accum.df<-rbind(accum.df,df)
  }
  
  return(accum.df)
}
