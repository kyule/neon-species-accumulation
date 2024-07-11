
### Want to understand years of sampling required, but sampling differs across years, need to scale by the number of traps

library(lubridate)
library(vegan)
library(ggplot2)
library(neonUtilities)


#### BEFORE RUNNING

#### Users must define their own paths and NEON token

# datapath<-"user defined path"
# codepath<-"user defined path"
# Neon_Token<-"user's token"
source("/Users/kelsey/Github/neon-species-accumulation/configini.R")

#### And users should define their variables for the NEON data download

# data product of interest, below is for Carabids
product<-"DP1.10022.001"

# start and end dates, must be formatted as "YYYY-MM" or NA for all time
start<-"2000-01"
end<-"2024-05"

# Use loadByProduct to download the data
NeonData<-loadByProduct(dpID=product,
                        site="PUUM",
                        startdate=start,
                        enddate=end,
                        token=Neon_Token, 
                        check.size=FALSE, 
                        include.provisional=FALSE)
save(NeonData, file=paste0("~/Desktop/","PUUM_NeonData.Robj"))

####
# Pull tables of interest
field<-NeonData$bet_fielddata
sort<-NeonData$bet_sorting
para<-NeonData$bet_parataxonomistID
expert<-NeonData$bet_expertTaxonomistIDProcessed

### Subset to only traps that were collected
field<-field[which(field$sampleCollected=="Y"),]

### Add year and summarize
field$year<-year(field$collectDate)
trapsPerYear<-field %>% group_by(year) %>% count()

randomSpAccum<-function(field.dat,traps.peryear,sort.dat,n.traps,n.iter){
  
  set.seed(85705)
  
  # subset to carabids
  sort.dat<-sort.dat[which(sort.dat$sampleType %in% c("carabid")),]
  
  # create data frame for species accumulation work
  accum.df<-data.frame(iter=c(),n.traps=c(),year=c(),richness=c())
  
  # cycle through the number of iterations
  for (i in 1:n.iter){
    
    #randomly reorder field data
    field.iter<-field[sample(1:nrow(field),nrow(field)),]
    
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
    rand.sort<-sort[which(sort$sampleID %in% rand.field$sampleID),]
    rand.sort$year<-year(rand.sort$collectDate)
    
    count.traps<-data.frame(rand.field %>% group_by(year) %>% count())$n
    
    rand.sort.totals<-data.frame(rand.sort %>% group_by(taxonID,year) %>% summarise(individualCount=sum(individualCount)))
    rand.sort.totals<-rand.sort.totals[!is.na(rand.sort.totals$individualCount),]
    #for (k in 1:nrow(trapsPerYear)){
      #year<-trapsPerYear$year[k]
      #rand.trap.no<-sum(rand.sort.totals$year==year)
      #if(rand.trap.no==0){rbind(rand.sort.tot)}
    #}

    
    comm<-reshape(rand.sort.totals[,c("taxonID","year","individualCount")],direction="wide",timevar="taxonID",idvar="year")
    comm[is.na(comm)]<-0

    print(comm)
    accum<-specaccum(comm)
    df<-data.frame(iter=rep(i,nrow(trapsPerYear)),n.traps=count.traps,year=accum$sites,richness=accum$richness)
    accum.df<-rbind(accum.df,df)
  }
  return(accum.df)
}

srer.25<-randomSpAccum(field,sort,25,100)
srer.100<-randomSpAccum(field,sort,100,100)
srer.max<-randomSpAccum(field,sort,min(trapsPerYear$n),100)
srer.obs<-randomSpAccum(field,sort,max(trapsPerYear$n),100)

ggplot(data=srer.obs[1:length(unique(srer.obs$year)),],aes(x=year,y=richness,group=iter,color='red'))+
  geom_line(linewidth=1.5)+
  theme_minimal()+
  theme(legend.position = 'none')+ 
  ylim(0, (max(srer.obs$richness)+10))+
  geom_line(data=srer.25,aes(x=year,y=richness,group=iter),color='lightgray',linewidth=0.25)+
  geom_line(data=srer.100,aes(x=year,y=richness,group=iter),color='gray',linewidth=0.25)+
  geom_line(data=srer.max,aes(x=year,y=richness,group=iter),color='darkgray',linewidth=0.25)

