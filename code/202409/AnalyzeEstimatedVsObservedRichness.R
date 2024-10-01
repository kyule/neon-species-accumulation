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


#### Start some basic analyses
summary(lm(Observed~turnover*years,full.com))
plot(Observed~turnover,full.com)
summary(lm(Observed~dissim*years,full.com))
plot(Observed~dissim,full.com)
# No relationship between observed and turnover or number of years / dissim number of years

summary(lm(Estimator~turnover*years,full.com))
plot(Estimator~turnover,full.com)
summary(lm(Estimator~dissim*years,full.com))
plot(Estimator~dissim,full.com)
# same for estimated


### Proportion observed analyses
ggplot(full.com, aes(x = turnover, y = propObs, color=as.numeric(Estimator))) +
  geom_point(aes(size=years)) + 
  geom_smooth(method = "glm", method.args = list(family = "quasibinomial"), color = "black") +  
  labs(x = "Mean Species Turnover", y = "Prop. Estimated Species Observed") +
  theme_minimal() +
  scale_color_viridis_c(option = "D",name="Est. Richness")
# negative relationship between turnover and proportion of estimated species richness we have observed


# Likelihood ratio test to determine best models of proportion observed
model_full <- glm(propObs ~ years * turnover * Estimator, family = quasibinomial, data = full.com)
model_reduced <- glm(propObs ~ years + turnover + Estimator + years:turnover + years:Estimator + turnover:Estimator, family = quasibinomial, data = full.com)
anova(model_reduced, model_full, test = "Chisq")
summary(model_reduced)

model_reduced2 <- glm(propObs ~ years + turnover + Estimator +  years:Estimator + turnover:Estimator, family = quasibinomial, data = full.com)
anova(model_reduced2, model_reduced, test = "Chisq")
summary(model_reduced2)

model_reduced3 <- glm(propObs ~ years + turnover + Estimator +  years:Estimator , family = quasibinomial, data = full.com)
anova(model_reduced3, model_reduced2, test = "Chisq")
summary(model_reduced3)

model_reduced4 <- glm(propObs ~ years + turnover + Estimator , family = quasibinomial, data = full.com)
anova(model_reduced4, model_reduced3, test = "Chisq")
summary(model_reduced3)
### The best model says that the proportion of estimated species that have been observed increases with number of years 
    #of sampling, decreases with average year-to-year species turnover, and 
    #decreases with the number of species that have been estimated

### Plot observed vs Estimated with errors around variance

# Determine overlap vs. no

full.com$signif<- ifelse(full.com$Observed > (full.com$Estimator + full.com$Est_s.e.) |
                           full.com$Observed < (full.com$Estimator - full.com$Est_s.e.), 
                         1, 0)

# Create rank order values
full.com<-full.com[order(full.com$Observed),]
full.com$obsRank<-1:nrow(full.com)
full.com<-full.com[order(full.com$turnover),]
full.com$turnRank<-1:nrow(full.com)
full.com<-full.com[order(full.com$dissim),]
full.com$dissimRank<-1:nrow(full.com)

ggplot(full.com, 
       aes(x = turnRank, y = Observed)) +
  geom_rect(aes(xmin = turnRank - 0.5, 
                xmax = turnRank + 0.5, 
                ymin = -Inf, 
                ymax = Inf, alpha = signif),
            fill = "grey") +
  geom_errorbar(aes(ymin = Estimator - Est_s.e., 
                    ymax = Estimator + Est_s.e., 
                    color = years), width = 0,size=2) + 
  geom_point(color = "black", size = 2) + 
  #geom_point(aes(x = 1:nrow(full.com), y = Estimator, color = propObs), size = 2) + 
  labs(x = "Rank-order Turnover", y = "Richness") +
  theme_bw() + 
  scale_color_viridis_c(option = "D", name = "years") +
  scale_alpha(range = c(0, 0.5), guide = "none")  

# Overlap with estimator +/- se
ggplot(full.com,aes(x=as.factor(signif),y=turnover))+
  geom_violin()


# Plot on "maps"
sites<-read.csv("/Users/kelsey/Github/neon-species-accumulation/data/NEON_Field_Site_Metadata_20240926.csv")
full.com<-left_join(full.com,sites,join_by("site"=="field_site_id"))

ggplot(full.com,aes(x=field_longitude,y=field_latitude)) +
  geom_point(aes(size=Estimator,color=as.numeric(propObs))) +
  theme_bw() + 
  scale_color_viridis_c(option = "D", name = "Obs/Exp")
  
# Adjusted pseudo-R2 model comparison

model <- glm(propObs ~ turnover + Estimator + years + Estimator:years, family = quasibinomial, data = full.com)
null_model <- glm(propObs ~ 1, family = quasibinomial, data = full.com)

dev_r2_mod <- 1 - (deviance(model) / deviance(null_model))
model2 <- glm(propObs ~ Estimator + years + Estimator:years , family = quasibinomial, data = full.com)
dev_r2_mod2 <- 1 - (deviance(model2) / deviance(null_model))

dev_r2_mod-dev_r2_mod2

n <- nrow(full.com)

# Number of predictors (excluding the intercept)
p1 <- length(coef(model)) - 1
p2 <- length(coef(model2)) - 1


# Adjusted pseudo-R-squared
adj_dev_r2_mod <- 1 - ((1 - dev_r2_mod) * (n - 1) / (n - p2 - 1))
adj_dev_r2_mod2 <- 1 - ((1 - dev_r2_mod2) * (n - 1) / (n - p2 - 1))
adj_dev_r2_mod - adj_dev_r2_mod2





