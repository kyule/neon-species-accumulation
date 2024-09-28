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
  source(paste0(codepath,"accumulation.R"))}else{load(file=paste0(datapath,"results.Robj"))}


# Remove GUAN due to very low sample size
results<- results[!names(results)=="GUAN"]

# Pull Estimated Asymptotic and Observed Richness Values out of the results list -- full data only

asyest.list<-lapply(results,function(item) data.frame(item$full$out$AsyEst))
asyest <- bind_rows(asyest.list,.id="site")
richness <- asyest[grep("Richness",rownames(asyest)),]

# Pull the turnovers out of the results list and input into the richness data frame

turnover.list<-lapply(results,function(item) data.frame(item$turnover))
turnover.list <- bind_rows(turnover.list,.id="site")
turnover<-data.frame(turnover.list %>% 
                       group_by(site) %>% 
                       summarise(turnover=mean(total,na.rm=TRUE),years=length(year)))
full.com <- left_join(richness,turnover,join_by("site"=="site"))


full.com$propObs<-full.com$Observed/full.com$Estimator

# Take into account the number of traps that have been sampled

trap.sum <- data.frame(completeness %>% group_by(siteID) %>% summarise(traps=sum(traps)))
full.com<-left_join(full.com,trap.sum,join_by("site"=="siteID"))

# Start some basic analyses
summary(lm(Observed~turnover,full.com))
plot(Observed~turnover,full.com)
# No relationship between observed and turnover

summary(lm(Observed~turnover,full.com))
plot(Observed~turnover,full.com)
# No relationship between observed and turnover

summary(lm(Estimator~turnover,full.com))
plot(Estimator~turnover,full.com)
# No relationship between estimated and turnover

summary(glm(propObs~turnover,full.com,family=binomial(link="logit")))

ggplot(full.com, aes(x = turnover, y = propObs, color=as.numeric(Estimator))) +
  geom_point(aes(size=traps)) + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), color = "black") +  # Add linear regression line
  labs(x = "Mean Species Turnover", y = "Prop. Estimated Species Observed") +
  theme_minimal() +
  scale_color_viridis_c(option = "D",name="Est. Richness")
# negative relationship between turnover and proportion of estimated species richness we have observed


summary(glm(propObs~Estimator,full.com,family=binomial(link="logit")))
# negative relationship between estimated spp richness and proportion of estimated species richness we have observed
ggplot(full.com, aes(x = Estimator, y = propObs, color=as.numeric(turnover))) +
  geom_point(size=3) + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), color = "black") +  # Add linear regression line
  labs(y = "Prop. Estimated Species Observed", x = "Estimated Richness") +
  theme_minimal() +
  scale_color_viridis_c(option = "D",name="turnover")
# negative relationship between the proportion observed and the estimated richness overall

summary(glm(propObs~years+turnover,full.com,family=binomial(link="logit")))
#There is no relationship between the proportion observed and the number of years sampled overall

summary(lm(Estimator~years,full.com))
# There is however a positive relationship between the estimated richness and the number of years of sampling
### Because turnover is high, estimated richness becomes larger the more years you sample

summary(glm(propObs~years*Estimator,full.com,family=binomial(link="logit")))
#Proportion observed decreases with the estimated number of species
# but that negative effect is decreased with increasing numbers of years of sampling

summary(lm(Estimator~traps,full.com))
plot(Estimator~traps,full.com)

### Find sampling required to reach 90% diversity

# Pull Estimated and Observed Richness by sampling Values out of the results list -- full community only

inext.list<-lapply(results,function(item) data.frame(item$full$out$iNextEst$size_based))
inext <- bind_rows(inext.list,.id="site")
inext <- inext %>% filter(Order.q==0)

thresh90<-data.frame(site=full.com$site,thresh=full.com$Estimator*0.90)
trapAvg<-data.frame(completeness %>% group_by(site) %>% )


ggplot(inext, aes(x = t, y = qD, color=site , group=site)) +
  geom_line() +
  labs(x = "Number of Traps", y = "Estimated Richness") +
  theme_minimal()






