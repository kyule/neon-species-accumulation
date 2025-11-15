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

load(file=paste0(datapath,"iNEXTandTurnoverResults_withSimpson.Robj"))
full.com<-read.csv(paste0(datapath,'communityResults_withSimpson.csv'))
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
inext <- inext %>% filter(Order.q==2)

thresh90<-data.frame(site=full.com$site,thresh=full.com$Estimator*0.90)
full.com$trapAvg<-full.com$traps/full.com$years
thresh90<-left_join(thresh90,full.com[,c("site","trapAvg")],join_by("site"=="site"))


## Find where the thresholds is crossed and scale years in inext data frame
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

# Bring in field site data, join other data for easy plotting

sites<-read.csv(paste0(datapath,"NEON_Field_Site_Metadata_20240926.csv"))
domains<-data.frame(ID=c('D01',	'D02',	'D03',	'D04',	'D05',	'D06',	'D07',	'D08',	'D09',	'D10',	'D11',	'D12',	'D13',	'D14',	'D15',	'D16',	'D17',	'D18',	'D19',	'D20'),
                         domainName=c('Northeast',	'Mid-Atlantic',	'Southeast',	'Atlantic Neotropical',	'Great Lakes',	'Prairie Peninsula',	'Central Plains',	'Ozarks Complex',	'Northern Plains',	'Southern Plains',	'Southern Rockies/Colorado Plateau',	'Desert Southwest',	'Pacific Southwest',	'Northern Rockies',	'Great Basin',	'Pacific Northwest',	'California',	'Tundra',	'Taiga',	'Pacific Tropical'))
sites<-left_join(sites,domains,join_by("field_domain_id"=="ID"))
inext<-left_join(inext,sites,join_by("site"=="field_site_id"))
thresh90<-left_join(thresh90,sites,join_by("site"=="field_site_id"))
thresh90<-left_join(thresh90,turnover,join_by("site"=="site"))

inext<-left_join(inext,thresh90,join_by("site"=="site"))


# Histogram colored by turnover

bin_stats <- thresh90 %>%
  mutate(bin = cut(y.thresh, breaks = c(0,1,2,3), right = TRUE)) %>%  
  group_by(bin) %>%
  summarise(
    count = n(),
    avg_z = mean(turnover, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(bin_center = (as.numeric(sub("\\((.+),.*", "\\1", bin)) + 
                         as.numeric(sub(".*,(.+)\\]", "\\1", bin))) / 2) 

turnover_limits <- range(c(thresh90$turnover), na.rm = TRUE)

bin <- ggplot(bin_stats, aes(x = bin, y = count, fill = avg_z)) +
  geom_bar(stat = "identity", width = 0.9) +  # smaller width so bars are not oversized
  scale_fill_viridis(option = "viridis", name = "Mean turnover", limits = turnover_limits) +
  labs(x = "Number of Years", y = "Number of sites") +
  theme_minimal()

print(bin)


# correlation with turnover

summary(thresh90$y.thresh)
sd(thresh90$y.thresh)/sqrt(nrow(thresh90))

cor.test(as.numeric(thresh90$y.thresh),as.numeric(thresh90$turnover))


##
div<-read.csv(paste0(datapath,'communityResults.csv'))

correlate<-left_join(full.com,div,join_by("site"=="site","year"=="year"))
cor(correlate$Observed.div,correlate$Observed)







