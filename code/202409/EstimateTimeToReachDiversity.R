### Find sampling required to reach 90% diversity

### This is the primary file for analyzing NEON species accumulation data

#### BEFORE RUNNING
# Users must define their own paths

# datapath<-"user defined path"
# codepath<-"user defined path"

# Users need to indicate whether they want to get new set of results
# These these steps are very time consuming so it is recommended that they are only done if necessary
NewResults<-FALSE

# And users should remove below line that loads in personal paths
source("/Users/kelsey/Github/neon-species-accumulation/configini.R")

#set seed
set.seed(85705)

# Load some necessary packages
library("dplyr")
library("ggplot2")
library("MASS")

# Load in the formatted clean data, or download and create it. 
#Make sure the results are correctly configured

if(NewResults==TRUE|file.exists(paste0(datapath,"results.Robj"))==FALSE){
  source(paste0(codepath,"accumulation.R"))}else{load(file=paste0(datapath,"resultsFull.Robj"))}

# Remove GUAN due to very low number of beetles captured. 
### Only 36 beetles in 7 years of sampling, almost complete turnover in communities between most years
results<- results[!names(results)=="GUAN"]

### Pull Estimated Asymptotic and Observed Richness Values out of the results list -- full data only

asyest.list<-lapply(results,function(item) data.frame(item$full$out$AsyEst))
asyest <- bind_rows(asyest.list,.id="site")
richness <- asyest[grep("Richness",rownames(asyest)),]

#### Pull the turnovers out of the results list and input into the richness data frame

turnover.list<-lapply(results,function(item) data.frame(item$turnover))
turnover.list <- bind_rows(turnover.list,.id="site")
turnover<-data.frame(turnover.list %>% 
                       group_by(site) %>% 
                       summarise(turnover=mean(total,na.rm=TRUE),years=length(year)))
full.com <- left_join(richness,turnover,join_by("site"=="site"))

full.com$propObs<-full.com$Observed/full.com$Estimator

# Take into account the number of traps that have been sampled
trap.list<-lapply(results,function(item) nrow(data.frame(item$full$traps)))
trap.list <- melt(trap.list)
names(trap.list)<-c("traps","site")
full.com<-left_join(full.com,trap.list,join_by("site"=="site"))

# Pull in dissimilarity measures
dissim<-melt(lapply(results,function(item) mean(item$dissim)))
names(dissim)<-c("dissim","site")
full.com<-left_join(full.com,dissim,join_by("site"=="site"))

# Pull Estimated and Observed Richness by sampling Values out of the results list -- full community only

inext.list<-lapply(results,function(item) data.frame(item$full$out$iNextEst$size_based))
inext <- bind_rows(inext.list,.id="site")
inext <- inext %>% filter(Order.q==0)

thresh90<-data.frame(site=full.com$site,thresh=full.com$Estimator*0.90)
full.com$trapAvg<-full.com$traps/full.com$years
thresh90<-left_join(thresh90,full.com[,c("site","trapAvg")],join_by("site"=="site"))

## Find where the threshold is crossed and scale years in inext data frame
thresh90$t.thresh<-NA
inext$y<-NA

for (i in 1:nrow(thresh90)){
  site_data<-inext[which(inext$site==thresh90$site[i]),]
  inext$y[which(inext$site==thresh90$site[i])]<-inext$t[which(inext$site==thresh90$site[i])]/thresh90$trapAvg[i]
  diff<-site_data$qD-thresh90$thresh[i]
  
  for (j in 1:length(diff)){
    if (diff[j]>0){break}
  }
  
  v<-site_data[c(j-1,j),c("t","qD")]
  thresh90$t.thresh[i]<-v$t[1] + (thresh90$thresh[i] - v$qD[1]) * (v$t[2] - v$t[1]) / (v$qD[2] - v$qD[1])
  
}


# Calculate the thresholds based on year instead 
thresh90$y.thresh<-thresh90$t.thresh/thresh90$trapAvg
full.com<-left_join(full.com,thresh90,join_by("site"=="site"))

# Bring in field site data

sites<-read.csv("/Users/kelsey/Github/neon-species-accumulation/data/NEON_Field_Site_Metadata_20240926.csv")






#Plot the accumulation curves and estimated number of years
ggplot(inext, aes(x = y, y = qD, color=site , group=site)) +
  geom_line() +
  labs(x = "Years", y = "Estimated Richness") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1") +
  geom_hline(data = thresh90, aes(yintercept = thresh, color = site), linetype = "dashed") +
  geom_point(data = thresh90, aes(x = y.thresh, y = thresh), color = "darkgrey", size = 2) 





# Scale the X axis by years of sampling

inext$byYears<-NA
for (i in 1:nrow(thresh90)){
  inext$byYears[which(inext$site==thresh90$site[i])]<-inext$t[which(inext$site==thresh90$site[i])]/thresh90$trapno[i]
}

intersections <- do.call(rbind, lapply(1:nrow(thresh90), function(i) {
  find_intersections(inext, thresh90[i, ],"byYears")
}))


ggplot(inext, aes(x = byYears, y = qD, color=site , group=site)) +
  geom_line() +
  labs(x = "Number of Years", y = "Estimated Richness") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1") +
  geom_hline(data = thresh90, aes(yintercept = thresh, color = site), linetype = "dashed") +
  geom_point(data = intersections, aes(x = t, y = qD), color = "black", size = 2) 



# Turnover relates to number of years?

thresh90<-left_join(thresh90,turnover,join_by("site"=="site"))

summary(lm(yearsReq~turnover,thresh90))






