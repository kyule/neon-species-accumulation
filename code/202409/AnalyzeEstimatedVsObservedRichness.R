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


# Remove GUAN due to very low 
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

trap.list<-lapply(results,function(item) nrow(data.frame(item$full$traps)))
trap.list <- melt(trap.list)
names(trap.list)<-c("traps","site")
full.com<-left_join(full.com,trap.list,join_by("site"=="site"))

# Start some basic analyses
summary(lm(Observed~turnover*years,full.com))
plot(Observed~turnover,full.com)
# No relationship between observed and turnover

summary(lm(Estimator~turnover,full.com))
plot(Estimator~turnover,full.com)
# Estimated diversity increases with species turnover


summary(glm(propObs~turnover,full.com,family=quasibinomial(link="logit")))

ggplot(full.com, aes(x = turnover, y = propObs, color=as.numeric(Estimator))) +
  geom_point(aes(size=traps)) + 
  geom_smooth(method = "glm", method.args = list(family = "quasibinomial"), color = "black") +  # Add linear regression line
  labs(x = "Mean Species Turnover", y = "Prop. Estimated Species Observed") +
  theme_minimal() +
  scale_color_viridis_c(option = "D",name="Est. Richness")
# negative relationship between turnover and proportion of estimated species richness we have observed


summary(glm(propObs~Estimator,full.com,family=quasibinomial(link="logit")))
# negative relationship between estimated spp richness and proportion of estimated species richness we have observed
ggplot(full.com, aes(x = Estimator, y = propObs, color=as.numeric(turnover))) +
  geom_point(size=3) + 
  geom_smooth(method = "glm", method.args = list(family = "quasibinomial"), color = "black") +  # Add linear regression line
  labs(y = "Prop. Estimated Species Observed", x = "Estimated Richness") +
  theme_minimal() +
  scale_color_viridis_c(option = "D",name="turnover")
# negative relationship between the proportion observed and the estimated richness overall

summary(glm(propObs~years+turnover,full.com,family=quasibinomial(link="logit")))


summary(lm(Estimator~years,full.com))
#estimated richness increases with years of sampling

fullmodel<-glm(propObs~years*Estimator*turnover*traps,full.com,family=quasibinomial(link="logit"))
summary(fullmodel)

library(caret)
train_control <- trainControl(method = "cv", number = 10)
model <- train(propObs~years*Estimator*turnover,full.com, method = "glmnet",family="quasibinomial", trControl = train_control)
summary(model)


# Fit two models: a simpler model and a more complex one
model_full <- glm(propObs ~ years * turnover * Estimator, family = quasibinomial, data = full.com)
model_reduced <- glm(propObs ~ years + turnover + Estimator + years:turnover + years:Estimator + turnover:Estimator, family = quasibinomial, data = full.com)
# Perform a likelihood ratio test
anova(model_reduced, model_full, test = "Chisq")
summary(model_reduced)

model_reduced2 <- glm(propObs ~ years + turnover + Estimator +  years:Estimator + turnover:Estimator, family = quasibinomial, data = full.com)
anova(model_reduced2, model_reduced, test = "Chisq")
summary(model_reduced2)

model_reduced3 <- glm(propObs ~ years + turnover + Estimator +  years:Estimator , family = quasibinomial, data = full.com)
anova(model_reduced3, model_reduced2, test = "Chisq")
summary(model_reduced3)

model_reduced4 <- glm(propObs ~ years + turnover + Estimator  , family = quasibinomial, data = full.com)
anova(model_reduced4, model_reduced3, test = "Chisq")
summary(model_reduced4)
### The best model says that the proportion of estimated species that have been observed increases with number of years of sampling, decreases with average year-to-year species turnover, and decreases with the number of species that have been estimated




summary(lm(Estimator~traps,full.com))
plot(Estimator~traps,full.com)
# The more traps you have total the more species you estimate


### Find sampling required to reach 90% diversity

# Pull Estimated and Observed Richness by sampling Values out of the results list -- full community only

inext.list<-lapply(results,function(item) data.frame(item$full$out$iNextEst$size_based))
inext <- bind_rows(inext.list,.id="site")
inext <- inext %>% filter(Order.q==0)

thresh90<-data.frame(site=full.com$site,thresh=full.com$Estimator*0.90)
full.com$trapAvg<-full.com$traps/full.com$years
thresh90<-left_join(thresh90,full.com[,c("site","trapAvg")],join_by("site"=="site"))

# Find intersections
find_intersections <- function(df,thresh,metric) {
  site_data <- df[df$site == thresh$site, ]
  
  # Loop through data and check for crossings
  for (i in 1:(nrow(site_data) - 1)) {
    if ((site_data$qD[i] < thresh$thresh && site_data$qD[i+1] > thresh$thresh) ||
        (site_data$qD[i] > thresh$thresh && site_data$qD[i+1] < thresh$thresh)) {
      
      # Linear interpolation to find the exact crossing point
      if (metric=="t"){
        t1 <- site_data$t[i]
        t2 <- site_data$t[i+1]
        qD1 <- site_data$qD[i]
        qD2 <- site_data$qD[i+1]
        crossing_t <- t1 + (thresh$thresh - qD1) * (t2 - t1) / (qD2 - qD1)
        
        return(data.frame(site = thresh$site, t = crossing_t, qD = thresh$thresh))
      }
      
      if (metric=="byYears")
      y1 <- site_data$byYears[i]
      y2 <- site_data$byYears[i+1]
      qD1 <- site_data$qD[i]
      qD2 <- site_data$qD[i+1]
      crossing_t <- y1 + (thresh$thresh - qD1) * (y2 - y1) / (qD2 - qD1)
      
      return(data.frame(site = thresh$site, t = crossing_t, qD = thresh$thresh))
    }
  }
  return(NULL)
}

# Get all intersections 
intersections <- do.call(rbind, lapply(1:nrow(thresh90), function(i) {
  find_intersections(inext, thresh90[i, ],"t")
}))




ggplot(inext, aes(x = t, y = qD, color=site , group=site)) +
  geom_line() +
  labs(x = "Number of Traps", y = "Estimated Richness") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1") +
  geom_hline(data = thresh90, aes(yintercept = thresh, color = site), linetype = "dashed") +
  geom_point(data = intersections, aes(x = t, y = qD), color = "darkgrey", size = 2) 

#Calculate the number of years based on recent per year sampling effort

thresh90<-left_join(thresh90,intersections[c("site","t")],join_by("site"=="site"))
thresh90$yearsReq<-thresh90$t/thresh90$trapno
thresh90

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






