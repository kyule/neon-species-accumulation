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
library("reshape")

# Load in the formatted clean data, or download and create it. 
#Make sure the results are correctly configured

#if(NewResults==TRUE|file.exists(paste0(datapath,"results.Robj"))==FALSE){source(paste0(codepath,"accumulation.R"))}else{load(file=paste0(datapath,"resultsFull.Robj"))}

# Remove GUAN due to very low number of beetles captured. 
### Only 36 beetles in 7 years of sampling, almost complete turnover in communities between most years
results<- results[!names(results)=="GUAN"]
#results<- results[!names(results)=="YELL"] take out in provisional
#results<- results[!names(results)=="TEAK"] take out in provisional

### Pull Estimated Asymptotic and Observed Richness Values out of the results list -- full data only

asyest.list<-lapply(results,function(item) data.frame(item$full$out$AsyEst))
asyest <- bind_rows(asyest.list,.id="site")
richness <- asyest[grep("Richness",rownames(asyest)),]

#### Pull the turnovers out of the results list and input into the richness data frame

turnover.list<-lapply(results,function(item) data.frame(item$turnover))
turnover.list <- bind_rows(turnover.list,.id="site")
turnover<-data.frame(turnover.list %>% 
                       group_by(site) %>% 
                       summarise(turnover=mean(total,na.rm=TRUE),years=length(year)+1))
full.com <- left_join(richness,turnover,join_by("site"=="site"))

full.com$propObs<-full.com$Observed/full.com$Estimator

# Take into account the number of traps that have been sampled
trap.list<-lapply(results,function(item) nrow(data.frame(item$full$traps)))
trap.list <- melt(trap.list)
names(trap.list)<-c("traps","site")
full.com<-left_join(full.com,trap.list,join_by("site"=="site"))

### Proportion observed analyses
ggplot(full.com, aes(x = turnover, y = propObs, color = as.numeric(Estimator))) +
  geom_point(aes(size = years)) + 
  geom_smooth(method = "glm", method.args = list(family = "quasibinomial"), color = "black") +  
  annotate("text", x = -Inf, y = Inf, label = "a", fontface = "bold", hjust = -25, vjust = 1.1, size = 6) +
  labs(x = "Mean Species Turnover", y = "Observed/Estimated Richness") +
  theme_minimal() +
  scale_color_viridis_c(option = "D", name = "Est. richness") +
  theme(
    axis.title = element_text(size = 18),  # Adjusts font size for axis titles
    axis.text = element_text(size = 18),   # Adjusts font size for axis labels
    legend.title = element_text(size = 18), # Adjusts font size for legend title
    legend.text = element_text(size = 18)   # Adjusts font size for legend text
  )
# negative relationship between turnover and proportion of estimated species richness we have observed


# Likelihood ratio test to determine best models of proportion observed
model_full <- glm(propObs ~ turnover *years *  Estimator, family = quasibinomial, data = full.com)
summary(model_full)
model_reduced <- glm(propObs ~ turnover + years + Estimator + turnover:years + turnover:Estimator + years:Estimator, family = quasibinomial, data = full.com)
anova(model_reduced, model_full, test = "Chisq")
summary(model_reduced)

model_reduced2 <- glm(propObs ~ turnover + years + Estimator + turnover:Estimator + years:Estimator, family = quasibinomial, data = full.com)
anova(model_reduced2, model_reduced, test = "Chisq")
summary(model_reduced2)

#with provisional and removing <5 years
model_reduced2 <- glm(propObs ~ turnover + years + Estimator + turnover:years + turnover:Estimator, family = quasibinomial, data = full.com)
model_reduced3 <- glm(propObs ~ turnover + years + Estimator + turnover:years, family = quasibinomial, data = full.com)
model_reduced4 <- glm(propObs ~ turnover + years + Estimator, family = quasibinomial, data = full.com)
model_reduced5 <- glm(propObs ~ turnover + Estimator, family = quasibinomial, data = full.com)
anova(model_reduced5, model_reduced4, test = "Chisq")
summary(model_reduced5)
model_reduced6 <- glm(propObs ~ Estimator, family = quasibinomial, data = full.com)
anova(model_reduced6, model_reduced5, test = "Chisq")
#Model 5 so decrease with estimator and turnover is the best

model_reduced3 <- glm(propObs ~ turnover + years + Estimator +  years:Estimator, family = quasibinomial, data = full.com)
anova(model_reduced3, model_reduced2, test = "Chisq")
summary(model_reduced3)

model_reduced4 <- glm(propObs ~ turnover + years + Estimator , family = quasibinomial, data = full.com)
anova(model_reduced4, model_reduced3, test = "Chisq")
summary(model_reduced4)

#stick with model 3
summary(model_reduced3)

#### The best model says that the proportion of estimated species that have been observed increases with number of years 
    #of sampling, decreases with average year-to-year species turnover, and 
    #decreases with the number of species that have been estimated

### Plot observed vs Estimated with errors around variance

# Determine overlap vs. no

full.com$signif<- ifelse(full.com$Observed > (full.com$Estimator + full.com$Est_s.e.) |
                           full.com$Observed < (full.com$Estimator - full.com$Est_s.e.), 
                         1, 0)

# Create rank order values
full.com<-full.com[order(full.com$Observed,decreasing=TRUE),]
full.com$obsRank<-1:nrow(full.com)
full.com<-full.com[order(full.com$turnover,decreasing=TRUE),]
full.com$turnRank<-1:nrow(full.com)


ggplot(full.com, 
       aes(x = obsRank, y = Observed)) +
  geom_rect(aes(xmin = obsRank - 0.5, 
                xmax = obsRank + 0.5, 
                ymin = -Inf, 
                ymax = Inf, alpha = signif),
            fill = "grey") +
  geom_errorbar(aes(ymin = Estimator - Est_s.e., 
                    ymax = Estimator + Est_s.e., 
                    color = turnover), width = 0,size=2) + 
  geom_point(color = "black", size = 2) + 
  #geom_point(aes(x = 1:nrow(full.com), y = Estimator, color = propObs), size = 2) + 
  labs(x = "Rank-order Observed Value", y = "Richness") +
  annotate("text", x = -Inf, y = Inf, label = "a", fontface = "bold", hjust = -30, vjust = 1.1, size = 6) +
  theme_bw() + 
  scale_color_viridis_c(option = "D", name = "turnover") +
  scale_alpha(range = c(0, 0.5), guide = "none")  +
  theme(
    axis.title = element_text(size = 18),  # Adjusts font size for axis titles
    axis.text = element_text(size = 18),   # Adjusts font size for axis labels
    legend.title = element_text(size = 18), # Adjusts font size for legend title
    legend.text = element_text(size = 18)   # Adjusts font size for legend text
  )

# Overlap with estimator +/- se
full.com$signifText<-"Overlap"
full.com$signifText[which(full.com$signif==1)]<-"No Overlap"
ggplot(full.com,aes(x=turnover, y=as.factor(signifText)))+
  labs(y = "Richness", x = "Mean Species Turnover") +
  geom_violin() + stat_summary(
    fun = "mean", geom = "point", shape = 20, size = 3, color = "black")+
  theme_minimal()   + theme(
    axis.title = element_text(size = 18),  # Adjusts font size for axis titles
    axis.text = element_text(size = 16),   # Adjusts font size for axis labels
  )
  

signif_model<-glm(signif~turnover,full.com,family="binomial")
summary(signif_model)

  