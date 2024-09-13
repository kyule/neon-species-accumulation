### This is the  file for downloading NEON data

#### BEFORE RUNNING

#### Users must define their own paths and NEON token

# datapath<-"user defined path"
# codepath<-"user defined path"
# Neon_Token<-"user's token"

#### And users should define their variables for the NEON data download



# load package
library(neonUtilities)


# Use loadByProduct to download the data
srer<-loadByProduct(dpID="DP1.10022.001",
                        site="SRER",
                        token=Neon_Token, 
                        check.size=FALSE, 
                        include.provisional=FALSE)

sort<-srer$bet_sorting
sort<-sort %>% filter(sampleType=="carabid")
field<-srer$bet_fielddata %>% filter(sampleCollected=="Y")



samps<-unique(field$sampleID)
spp<-unique(sort$taxonID)

inc<-data.frame(matrix(ncol=length(samps),nrow=length(spp)))
colnames(inc)<-samps
rownames(inc)<-spp

for (i in 1:nrow(inc)){
  for (j in 1:ncol(inc)){
    inds<-sort$individualCount[which(sort$sampleID==samps[j] & sort$taxonID==spp[i])]
    if (length(inds)>0) {inc[i,j]<-sum(inds)}
  }
}

inc[is.na(inc)]<-0
rowSums(inc)
colSums(inc)

presabs<-inc
presabs[presabs>1]<-1

# Convert to Incidence_freq
# incidence_freq
# freq <- c(our sample size, sp1= saw in 3 samples, sp2=saw in 8 samples)
# freq <- c(8, 3, 8)
input<-c(ncol(presabs),as.vector(rowSums(presabs)))

# input into plotsamps
out<-iNEXT(input,datatype='incidence_freq',knot=10)
#Because many traps don't contain any carabids we can't use it

out$DataInfo
out$AsyEst
out$iNextEst


###
sort<-left_join(sort,field,join_by(sampleID==sampleID))

events<-unique(sort$eventID)
spp<-unique(sort$taxonID)

inc<-data.frame(matrix(ncol=length(events),nrow=length(spp)))
colnames(inc)<-events
rownames(inc)<-spp

for (i in 1:nrow(inc)){
  for (j in 1:ncol(inc)){
    inds<-sort$individualCount[which(sort$eventID==events[j] & sort$taxonID==spp[i])]
    if (length(inds)>0) {inc[i,j]<-sum(inds)}
  }
}

inc[is.na(inc)]<-0
presabs<-inc
presabs[presabs>1]<-1

# Convert to Incidence_freq
input<-c(ncol(presabs),as.vector(rowSums(presabs)))

# input into inext
out<-iNEXT(input,datatype='incidence_freq',q=c(0,1,2),knot=10)

out$DataInfo
out$iNextEst
out$AsyEst

ggiNEXT(out)
ggiNEXT(out,type=2)
ggiNEXT(out,facet.var="Order.q",type=3)

### think about how granular to be given unequal trapping effort within events across years
### think about separating by year
