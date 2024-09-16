# This is the code for calculating dissimilarity estimates across years

library("rarestR")
library("corrplot")
library("reshape2")

#### BEFORE RUNNING
# Users must define their own paths

# datapath<-"user defined path"
# codepath<-"user defined path"

# And users must have configured the DownloadNEONData.R file as desired for their work
# Indicate whether new data should be downloaded, downloads are time consuming so this is not recommended
NewCleanData<-FALSE

# And users should remove below line that loads in personal paths
source("/Users/kelsey/Github/neon-species-accumulation/configini.R")

#set seed
set.seed(85705)

# Load in the formatted clean data, or download and create it. 
#Make sure DataCleaning.R is correctly configured

if(NewCleanData==TRUE|file.exists(paste0(datapath,"202409_CleanedData.Robj"))==FALSE){
  source(paste0(codepath,"DataCleaning.R"))}else{load(file=paste0(datapath,"202409_CleanedData.Robj"))}

# Define data of interest

data<-CleanedData$fullData

# Create community matrix for Year X Site 

spp<-unique(data$taxonID)
data$siteyear<-paste(data$domainID,data$siteID,data$year,sep="_")
siteyear<-unique(data$siteyear)
data<-data[order(data$siteyear),]

input<-data.frame(matrix(nrow=length(siteyear),ncol=length(spp)))
rownames(input)<-siteyear
colnames(input)<-spp

for (i in 1:nrow(input)){
  for (j in 1:ncol(input)){
    inds<-data$individualCountFinal[which(data$siteyear==rownames(input)[i]
                     & data$taxonID==colnames(input)[j])]
    if(length(inds)>0){input[i,j]<-sum(inds)}
    else(input[i,j]<-0)
  }
}

# Calculate and plot dissimilarities
dissim <- as.matrix(ess(input))
corrplot(dissim,method='color',col.lim=c(0,1))

# For a specific group of samples
factor<-"D16"
corDis<-dissim[grep(factor,rownames(input)),grep(factor,rownames(input))]
corrplot(corDis,method='color',col.lim=c(0,1))

# Turn data frame into long form
dissim.long<-melt(dissim)
names(dissim.long)<-c('samp1','samp2','dissim')
dissim.long$dom1<-substr(dissim.long$samp1,1,3)
dissim.long$dom2<-substr(dissim.long$samp2,1,3)
dissim.long$site1<-substr(dissim.long$samp1,5,8)
dissim.long$site2<-substr(dissim.long$samp2,5,8)
dissim.long$year1<-as.numeric(substr(dissim.long$samp1,10,13))
dissim.long$year2<-as.numeric(substr(dissim.long$samp2,10,13))
dissim.long$sameDomain<-0
dissim.long$sameDomain[which(dissim.long$dom1==dissim.long$dom2)]<-1
dissim.long$sameSite<-0
dissim.long$sameSite[which(dissim.long$site1==dissim.long$site2)]<-1
dissim.long$years<-abs(dissim.long$year1-dissim.long$year2)

# Plot Values

plot(dissim~as.factor(sameSite),dissim.long)
plot(dissim~as.factor(sameDomain),dissim.long)
plot(dissim~as.factor(years),dissim.long[which(dissim.long$sameSite==1),])









