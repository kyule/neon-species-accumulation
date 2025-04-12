### Find sampling required to reach 90% diversity

### This is the primary file for analyzing NEON species accumulation data

#### BEFORE RUNNING
# Users must define their own paths

# datapath<-"user defined path"
# codepath<-"user defined path"

# And users should remove below line that loads in personal paths
source("/Users/kelsey/Github/neon-species-accumulation/configini.R")

#set seed
set.seed(85705)

# Load some necessary packages
library("dplyr")
library("reshape")
library("ggplot2")
library("patchwork")
library('viridis')

# Load in the formatted clean data
#Make sure the results are correctly configured

load(file=paste0(datapath,"iNEXTandTurnoverResults.Robj"))
full.com<-read.csv(paste0(datapath,'communityResults.csv'))
full.com<-full.com[which(full.com$year=="full"),]

# Remove GUAN due to low data availability
results<- results[!names(results)=="GUAN"]

# Grab turnover

turnover.list<-lapply(results,function(item) data.frame(item$turnover))
turnover.list <- bind_rows(turnover.list,.id="site")
turnover<-data.frame(turnover.list %>% 
                       group_by(site) %>% 
                       summarise(turnover=mean(total,na.rm=TRUE),years=length(year)))


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
thresh90.rich<-left_join(thresh90.rich,turnover,join_by("site"=="site"))
thresh90.div<-left_join(thresh90.div,turnover,join_by("site"=="site"))

inext.rich<-left_join(inext.rich,thresh90.rich,join_by("site"=="site"))
inext.div<-left_join(inext.div,thresh90.div,join_by("site"=="site"))


# Histogram colored by turnover

bin_stats_rich <- thresh90.rich %>%
  mutate(bin = cut(y.thresh, breaks = 12)) %>%  
  group_by(bin) %>%
  summarise(count = n(), avg_z = mean(turnover, na.rm = TRUE), .groups = "drop") %>%
  mutate(bin_center = (as.numeric(sub("\\((.+),.*", "\\1", bin)) + 
                         as.numeric(sub(".*,(.+)\\]", "\\1", bin))) / 2) 

bin_stats_div <- thresh90.div %>%
  mutate(bin = cut(y.thresh, breaks = 12)) %>%  
  group_by(bin) %>%
  summarise(count = n(), avg_z = mean(turnover, na.rm = TRUE), .groups = "drop") %>%
  mutate(bin_center = (as.numeric(sub("\\((.+),.*", "\\1", bin)) + 
                         as.numeric(sub(".*,(.+)\\]", "\\1", bin))) / 2) 


turnover_limits<-range(c(thresh90.rich$turnover,thresh90.div$turnover),na.rm=TRUE)

rich_bin <- ggplot(bin_stats_rich, aes(x = bin_center, y = count, fill = avg_z)) +
  geom_bar(stat = "identity", width = 4) +
  scale_fill_viridis(option = "viridis", name = "Mean turnover", limits = turnover_limits) +
  labs(x = "", y = "Number of sites") +
  #annotate("text", x = -Inf, y = Inf, label = "b", fontface = "bold", hjust = -0.2, vjust = 1.3, size = 6) +
  theme_minimal()

div_bin <- ggplot(bin_stats_div, aes(x = bin_center, y = count, fill = avg_z)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_fill_viridis(option = "viridis", name =  "Mean turnover", limits = turnover_limits) +
  labs(x = "Years to reach threshhold", y = "Number of sites") +
  #annotate("text", x = -Inf, y = Inf, label = "d", fontface = "bold", hjust = -0.2, vjust = 1.3, size = 6) +
  theme_minimal()

combined_plot <- rich_bin / div_bin + 
  plot_layout(guides = "collect") + 
  theme(
    plot.caption = element_text(size = 18, hjust = 0.5, vjust = 1) 
  )

print(combined_plot)

# correlation with turnover

summary(rich.thresh$y.thresh)
sd(rich.thresh$y.thresh)/sqrt(nrow(rich.thresh))
summary(div.thresh$y.thresh)
sd(div.thresh$y.thresh)/sqrt(nrow(div.thresh))

cor.test(rich.thresh$y.thresh,rich.thresh$turnover)
cor.test(div.thresh$y.thresh,div.thresh$turnover)

# save threshold results

write.csv(thresh90.rich,paste0(datapath,'richnessThresh.csv'),row.names=FALSE)
write.csv(thresh90.div,paste0(datapath,'diversityThresh.csv'),row.names=FALSE)













