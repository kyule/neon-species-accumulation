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

# Load in the formatted clean data, or download and create it. 
#Make sure the results are correctly configured

if(NewResults==TRUE|file.exists(paste0(datapath,"results.Robj"))==FALSE){
  source(paste0(codepath,"accumulation.R"))}else{load(file=paste0(datapath,"results.Robj"))}

# Pull Estimated Asymptotic and Observed Richness Values out of the results list

asyest.list<-lapply(results,function(item) data.frame(item$out$AsyEst))
asyest <- bind_rows(asyest.list,.id="site")
richness <- asyest %>% filter(Diversity=="Species richness")
richness$Assemblage[which(richness$Assemblage=="full")]<-"0"
richness$Assemblage<-as.numeric(richness$Assemblage)

# Pull the turnovers out of the results list and input into the richness data frame

turnover.list<-lapply(results,function(item) data.frame(item$turnover))
turnover <- bind_rows(turnover.list,.id="site")
richness <- left_join(richness,turnover,join_by("site"=="site","Assemblage"=="year"))
richness<- richness %>% dplyr::rename(year=Assemblage) %>% dplyr::rename(turnover=total)

# Full community only analysis: Richness Vs Turnover
full.com<-richness[which(richness$year==0),]
# Give Average +/- se turnover per site
mean.turns<-richness %>% group_by(site) %>% 
  summarise(turn=mean(turnover,na.rm=TRUE),
            sd.turn = ifelse(sum(!is.na(turnover)) > 1, sd(turnover, na.rm = TRUE), NA)
  )

full.com$sd.turn<-NA

for (i in 1: nrow(mean.turns)){
  full.com$turnover[which(full.com$year==0 & full.com$site==mean.turns$site[i])]<-mean.turns$turn[i]
  full.com$sd.turn[which(full.com$year==0 & full.com$site==mean.turns$site[i])]<-mean.turns$sd.turn[i]
  
}

summary(lm(Observed~turnover,full.com))
plot(Observed~turnover,full.com)
# No relationship between observed and turnover

summary(lm(Estimator~turnover,full.com))
plot(Estimator~turnover,full.com)
# No relationship between estimated and turnover

hist(full.com$Observed/full.com$Estimator)
summary(lm(Observed/Estimator~turnover,full.com))
ggplot(full.com, aes(x = turnover, y = Observed / Estimator, color=as.numeric(Estimator))) +
  geom_point(size=3) + 
  geom_smooth(method = "lm", se = FALSE, color = "black") +  # Add linear regression line
  labs(x = "Mean Species Turnover", y = "Prop. Estimated Species Observed") +
  theme_minimal() +
  scale_color_viridis_c(option = "D",name="Est. Richness")
# negative relationship between turnover and proportion of estimated species richness we have observed


summary(lm(Observed/Estimator~Estimator,full.com))

ggplot(full.com, aes(x = Estimator, y = Observed/Estimator, color=as.numeric(turnover))) +
  geom_point(size=3) + 
  geom_smooth(method = "lm", se = FALSE, color = "black") +  # Add linear regression line
  labs(y = "Prop. Estimated Species Observed", x = "Estimated Richness") +
  theme_minimal() +
  scale_color_viridis_c(option = "D",name="turnover")
# Positive relationship between observed richness to date, and richness expected to not have been observed yet



