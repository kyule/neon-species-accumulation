### Find sampling required to reach 95% diversity

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

#if(NewResults==TRUE|file.exists(paste0(datapath,"results.Robj"))==FALSE){source(paste0(codepath,"accumulation.R"))}else{load(file=paste0(datapath,"resultsFull.Robj"))}

# Remove GUAN due to very low number of beetles captured. 
### Only 36 beetles in 7 years of sampling, almost complete turnover in communities between most years
results<- results[!names(results)=="GUAN"]

### Pull Estimated Asymptotic and Observed Diversity Values out of the results list -- full data only

asyest.list<-lapply(results,function(item) data.frame(item$full$out$AsyEst))
asyest <- bind_rows(asyest.list,.id="site")
diversity <- asyest[grep("Shannon",rownames(asyest)),]

#### Pull the turnovers out of the results list and input into the richness data frame

turnover.list<-lapply(results,function(item) data.frame(item$turnover))
turnover.list <- bind_rows(turnover.list,.id="site")
turnover<-data.frame(turnover.list %>% 
                       group_by(site) %>% 
                       summarise(turnover=mean(total,na.rm=TRUE),years=length(year)))
full.com <- left_join(diversity,turnover,join_by("site"=="site"))

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

# Pull Estimated and Observed Diversity by sampling Values out of the results list -- full community only

inext.list<-lapply(results,function(item) data.frame(item$full$out$iNextEst$size_based))
inext <- bind_rows(inext.list,.id="site")
inext <- inext %>% filter(Order.q==0)

thresh95<-data.frame(site=full.com$site,thresh=full.com$Estimator*0.95)
full.com$trapAvg<-full.com$traps/full.com$years
thresh95<-left_join(thresh95,full.com[,c("site","trapAvg")],join_by("site"=="site"))

## Find where the threshold is crossed and scale years in inext data frame
thresh95$t.thresh<-NA
inext$y<-NA

for (i in 1:nrow(thresh95)){
  site_data<-inext[which(inext$site==thresh95$site[i]),]
  inext$y[which(inext$site==thresh95$site[i])]<-inext$t[which(inext$site==thresh95$site[i])]/thresh95$trapAvg[i]
  diff<-site_data$qD-thresh95$thresh[i]
  
  for (j in 1:length(diff)){
    if (diff[j]>0){break}
  }
  
  v<-site_data[c(j-1,j),c("t","qD")]
  thresh95$t.thresh[i]<-v$t[1] + (thresh95$thresh[i] - v$qD[1]) * (v$t[2] - v$t[1]) / (v$qD[2] - v$qD[1])
  
}


# Calculate the thresholds based on year instead 
thresh95$y.thresh<-thresh95$t.thresh/thresh95$trapAvg
full.com<-left_join(full.com,thresh95,join_by("site"=="site"))

# Bring in field site data, join other data for easy plotting

sites<-read.csv("/Users/kelsey/Github/neon-species-accumulation/data/NEON_Field_Site_Metadata_20240926.csv")
domains<-data.frame(ID=c('D01',	'D02',	'D03',	'D04',	'D05',	'D06',	'D07',	'D08',	'D09',	'D10',	'D11',	'D12',	'D13',	'D14',	'D15',	'D16',	'D17',	'D18',	'D19',	'D20'),
                    domainName=c('Northeast',	'Mid-Atlantic',	'Southeast',	'Atlantic Neotropical',	'Great Lakes',	'Prairie Peninsula',	'Central Plains',	'Ozarks Complex',	'Northern Plains',	'Southern Plains',	'Southern Rockies/Colorado Plateau',	'Desert Southwest',	'Pacific Southwest',	'Northern Rockies',	'Great Basin',	'Pacific Northwest',	'California',	'Tundra',	'Taiga',	'Pacific Tropical'))
sites<-left_join(sites,domains,join_by("field_domain_id"=="ID"))
inext<-left_join(inext,sites,join_by("site"=="field_site_id"))
thresh95<-left_join(thresh95,sites,join_by("site"=="field_site_id"))
inext<-left_join(inext,full.com,join_by("site"=="site"))

# Find x and ylims
xlim<-ceiling(max(thresh95$y.thresh))
ylim<-ceiling(max(asyest$Estimator))

#Plot the accumulation curves and estimated number of years
ggplot(inext, aes(x = y, y = qD,group=site,color=as.numeric(turnover))) +
  geom_line(size=1) +
  facet_wrap(~ domainName) +
  labs(x = "Years", y = "Estimated Diversity") +
  theme_minimal() +
  geom_hline(data = thresh95, aes(yintercept = thresh), linetype = "dashed") +
  geom_point(data = thresh95, aes(x = y.thresh, y = thresh), color = "darkgrey", size = 2) +
  geom_point(data = inext[which(inext$Method=="Observed"),],aes(x= y, y=qD),color="black", size = 2) +
  ylim(0,ylim) +
  xlim(0,xlim) +
  scale_color_viridis_c(option = "D",name="Avg. turnover")

hist(thresh95$y.thresh)
summary(thresh95$y.thresh)

ggplot(thresh95, aes(x = y.thresh)) +
  geom_histogram(binwidth = 5,color='black') +
  scale_fill_viridis_c(option = "B") +
  labs(x = "Years of sampling to reach threshold", y = "Count") +
  theme_minimal()










