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
library("patchwork")

# Load in the formatted clean data, or download and create it. 
#Make sure the results are correctly configured

load(file=paste0(datapath,"iNEXTandTurnoverResults.Robj"))

# Remove GUAN due to very low number of beetles captured. 
### Only 36 beetles in 7 years of sampling, almost complete turnover in communities between most years
results<- results[!names(results)=="GUAN"]

### Pull Estimated Asymptotic and Observed Richness and diversity Values out of the results list -- full data only

asyest.list<-lapply(results,function(item) data.frame(item$full$out$AsyEst))
asyest <- bind_rows(asyest.list,.id="site")
richness <- asyest[grep("Richness",rownames(asyest)),]
diversity <- asyest[grep("Shannon",rownames(asyest)),]

full.com <- left_join(richness,diversity,join_by("site"=="site"),suffix = c(".rich",".div"))


#### Pull the turnovers out of the results list and input into the richness data frame

turnover.list<-lapply(results,function(item) data.frame(item$turnover))
turnover.list <- bind_rows(turnover.list,.id="site")
turnover<-data.frame(turnover.list %>% 
                       group_by(site) %>% 
                       summarise(turnover=mean(total,na.rm=TRUE),years=length(year)))
full.com <- left_join(full.com,turnover,join_by("site"=="site"))

full.com$propObs.rich<-full.com$Observed.rich/full.com$Estimator.rich
full.com$propObs.div<-full.com$Observed.div/full.com$Estimator.div

# Take into account the number of traps that have been sampled
trap.list<-lapply(results,function(item) nrow(data.frame(item$full$traps)))
trap.list <- melt(trap.list)
names(trap.list)<-c("traps","site")
full.com<-left_join(full.com,trap.list,join_by("site"=="site"))

# Pull Estimated and Observed Richness and diversity by sampling Values out of the results list -- full community only

inext.list<-lapply(results,function(item) data.frame(item$full$out$iNextEst$size_based))
inext <- bind_rows(inext.list,.id="site")
inext.rich <- inext %>% filter(Order.q==0)
inext.div <- inext %>% filter(Order.q==1)

thresh90.rich<-data.frame(site=full.com$site,thresh=full.com$Estimator.rich*0.90)
thresh90.div<-data.frame(site=full.com$site,thresh=full.com$Estimator.div*0.90)
full.com$trapAvg<-full.com$traps/full.com$years
thresh90.rich<-left_join(thresh90.rich,full.com[,c("site","trapAvg")],join_by("site"=="site"))
thresh90.div<-left_join(thresh90.div,full.com[,c("site","trapAvg")],join_by("site"=="site"))


## Find where the thresholds is crossed and scale years in inext data frame
thresh90.rich$t.thresh<-NA
inext.rich$y<-NA

for (i in 1:nrow(thresh90.rich)){
  site_data<-inext.rich[which(inext.rich$site==thresh90.rich$site[i]),]
  inext.rich$y[which(inext.rich$site==thresh90.rich$site[i])]<-inext.rich$t[which(inext.rich$site==thresh90.rich$site[i])]/thresh90.rich$trapAvg[i]
  diff<-site_data$qD-thresh90.rich$thresh[i]
  
  for (j in 1:length(diff)){
    if (diff[j]>0){break}
  }
  
  v<-site_data[c(j-1,j),c("t","qD")]
  thresh90.rich$t.thresh[i]<-v$t[1] + (thresh90.rich$thresh[i] - v$qD[1]) * (v$t[2] - v$t[1]) / (v$qD[2] - v$qD[1])
  
}

thresh90.div$t.thresh<-NA
inext.div$y<-NA

for (i in 1:nrow(thresh90.div)){
  site_data<-inext.div[which(inext.div$site==thresh90.div$site[i]),]
  inext.div$y[which(inext.div$site==thresh90.div$site[i])]<-inext.div$t[which(inext.div$site==thresh90.div$site[i])]/thresh90.div$trapAvg[i]
  diff<-site_data$qD-thresh90.div$thresh[i]
  
  for (j in 1:length(diff)){
    if (diff[j]>0){break}
  }
  
  v<-site_data[c(j-1,j),c("t","qD")]
  thresh90.div$t.thresh[i]<-v$t[1] + (thresh90.div$thresh[i] - v$qD[1]) * (v$t[2] - v$t[1]) / (v$qD[2] - v$qD[1])
  
}


# Calculate the thresholds based on year instead 
thresh90.rich$y.thresh<-thresh90.rich$t.thresh/thresh90.rich$trapAvg
thresh90.div$y.thresh<-thresh90.div$t.thresh/thresh90.div$trapAvg

full.com.rich<-left_join(full.com,thresh90.rich,join_by("site"=="site"))
full.com.div<-left_join(full.com,thresh90.div,join_by("site"=="site"))

# Bring in field site data, join other data for easy plotting

sites<-read.csv(paste0(datapath,"NEON_Field_Site_Metadata_20240926.csv"))
domains<-data.frame(ID=c('D01',	'D02',	'D03',	'D04',	'D05',	'D06',	'D07',	'D08',	'D09',	'D10',	'D11',	'D12',	'D13',	'D14',	'D15',	'D16',	'D17',	'D18',	'D19',	'D20'),
                         domainName=c('Northeast',	'Mid-Atlantic',	'Southeast',	'Atlantic Neotropical',	'Great Lakes',	'Prairie Peninsula',	'Central Plains',	'Ozarks Complex',	'Northern Plains',	'Southern Plains',	'Southern Rockies/Colorado Plateau',	'Desert Southwest',	'Pacific Southwest',	'Northern Rockies',	'Great Basin',	'Pacific Northwest',	'California',	'Tundra',	'Taiga',	'Pacific Tropical'))
sites<-left_join(sites,domains,join_by("field_domain_id"=="ID"))
inext<-left_join(inext,sites,join_by("site"=="field_site_id"))
thresh90.rich<-left_join(thresh90.rich,sites,join_by("site"=="field_site_id"))
thresh90.div<-left_join(thresh90.div,sites,join_by("site"=="field_site_id"))

inext.rich<-left_join(inext.rich,thresh90.rich,join_by("site"=="site"))
inext.div<-left_join(inext.div,thresh90.div,join_by("site"=="site"))

inext.rich<-left_join(inext.rich,turnover,join_by("site"=="site"))
inext.div<-left_join(inext.div,turnover,join_by("site"=="site"))

# Find x and ylims
xlim.rich<-ceiling(max(thresh90.rich$y.thresh))
ylim.rich<-ceiling(max(richness$Estimator))
xlim.div<-ceiling(max(thresh90.div$y.thresh))
ylim.div<-ceiling(max(diversity$Estimator))

#Plot the accumulation curves and estimated number of years
ggplot(inext.rich, aes(x = y, y = qD,group=site,color=as.numeric(turnover))) +
  geom_line(size=1) +
  facet_wrap(~ domainName) +
  labs(x = "Years", y = "Estimated Richness") +
  theme_minimal() +
  geom_hline(data = thresh90.rich, aes(yintercept = thresh), linetype = "dashed") +
  geom_point(data = thresh90.rich, aes(x = y.thresh, y = thresh), color = "darkgrey", size = 2) +
  geom_point(data = inext.rich[which(inext.rich$Method=="Observed"),],aes(x= y, y=qD),color="black", size = 2) +
  ylim(0,ylim.rich) +
  xlim(0,xlim.rich) +
  scale_color_viridis_c(option = "D",name="Avg. turnover")

ggplot(inext.div, aes(x = y, y = qD,group=site,color=as.numeric(turnover))) +
  geom_line(size=1) +
  facet_wrap(~ domainName) +
  labs(x = "Years", y = "Estimated Diversity") +
  theme_minimal() +
  geom_hline(data = thresh90.div, aes(yintercept = thresh), linetype = "dashed") +
  geom_point(data = thresh90.div, aes(x = y.thresh, y = thresh), color = "darkgrey", size = 2) +
  geom_point(data = inext.div[which(inext.div$Method=="Observed"),],aes(x= y, y=qD),color="black", size = 2) +
  ylim(0,ylim.div) +
  xlim(0,xlim.div) +
  scale_color_viridis_c(option = "D",name="Avg. turnover")


# Basic histogram
plotrich <- ggplot(thresh90.rich, aes(x = y.thresh)) +
  geom_histogram(binwidth = 5, color = 'black') +
  scale_y_continuous(limits = c(0, 20)) + 
  labs(y="Number of sites")+
  theme_minimal() +
  theme(
    axis.title.y = element_text(size = 18),
    axis.text = element_text(size = 18),
    axis.title.x = element_blank(),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18)
  )

plotdiv<-ggplot(thresh90.div, aes(x = y.thresh)) +
  geom_histogram(binwidth = 0.5,color='black') +
  scale_y_continuous(limits = c(0, 20)) +
  scale_x_continuous(limits = c(0, 8)) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_blank(),   
    axis.title.y = element_blank(),
    legend.title = element_text(size = 18), 
    legend.text = element_text(size = 18)   
  )

combined_plot <- plotrich + plotdiv + 
  plot_layout(guides = "collect") + 
  plot_annotation(
    caption = "Years to reach threshold"  
  ) & 
  theme(
    plot.caption = element_text(size = 18, hjust = 0.5, vjust = 1)  # Center and adjust position of caption
  )

print(combined_plot)










