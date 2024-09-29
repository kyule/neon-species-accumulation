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

### Pull Estimated Asymptotic and Observed Diversity Values out of the results list -- full data only

asyest.list<-lapply(results,function(item) data.frame(item$full$out$AsyEst))
asyest <- bind_rows(asyest.list,.id="site")
Diversity <- asyest[grep("Shannon",rownames(asyest)),]

#### Pull the turnovers out of the results list and input into the Diversity data frame

turnover.list<-lapply(results,function(item) data.frame(item$turnover))
turnover.list <- bind_rows(turnover.list,.id="site")
turnover<-data.frame(turnover.list %>% 
                       group_by(site) %>% 
                       summarise(turnover=mean(total,na.rm=TRUE),years=length(year)))
full.com <- left_join(Diversity,turnover,join_by("site"=="site"))

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
  labs(x = "Mean Species Turnover", y = "Q=1: obs/est") +
  theme_minimal() +
  scale_color_viridis_c(option = "D",name="Est. Diversity")+
  ylim(0.25,1)
# negative relationship between turnover and proportion of estimated species Diversity we have observed


# Likelihood ratio test to determine best models of proportion observed
model_full <- glm(propObs ~ years * turnover * Estimator, family = quasibinomial, data = full.com)
model_reduced <- glm(propObs ~ years + turnover + Estimator + years:turnover + years:Estimator + turnover:Estimator, family = quasibinomial, data = full.com)
anova(model_reduced, model_full, test = "Chisq")
summary(model_full)
### The best model includes all factors and all interactions: turnover has strong neg affect, year and estimator also have neg effects

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
  labs(x = "Rank-order Turnover", y = "Diversity") +
  theme_bw() + 
  scale_color_viridis_c(option = "D", name = "years") +
  scale_alpha(range = c(0, 0.5), guide = "none")  




