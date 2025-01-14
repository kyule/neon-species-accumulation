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
library("ggplot2")
library("MASS")
library("reshape")
library("cowplot")
library("gridExtra")
library("patchwork")

# Load in the formatted clean data, or download and create it. 
#Make sure the results are correctly configured

load(file=paste0(datapath,"iNEXTandTurnoverResults.Robj"))

# Remove GUAN due to very low number of beetles captured. 
### Only 36 beetles in 7 years of sampling, almost complete turnover in communities between most years
results<- results[!names(results)=="GUAN"]

### Pull Estimated Asymptotic and Observed Richness and Diversity Values out of the results list -- full data only

asyest.list<-lapply(results,function(item) data.frame(item$full$out$AsyEst))
asyest <- bind_rows(asyest.list,.id="site")
Diversity <- asyest[grep("Shannon",rownames(asyest)),]
Diversity<- Diversity[,1:4]
Richness <- asyest[grep("Richness",rownames(asyest)),]
Richness <- Richness[,1:4]

full.com <- left_join(Richness,Diversity,join_by("site"=="site"),suffix = c(".rich",".div"))

#### Pull the turnovers out of the results list and input into the full community data frame

turnover.list<-lapply(results,function(item) data.frame(item$turnover))
turnover.list <- bind_rows(turnover.list,.id="site")
turnover<-data.frame(turnover.list %>% 
                       group_by(site) %>% 
                       summarise(turnover=mean(total,na.rm=TRUE),years=length(year)+1))
full.com <- left_join(full.com,turnover,join_by("site"=="site"))

full.com$propObs.rich<-full.com$Observed.rich/full.com$Estimator.rich
full.com$propObs.div<-full.com$Observed.div/full.com$Estimator.div

# Take into account the number of traps that have been sampled
trap.list<-lapply(results,function(item) nrow(data.frame(item$full$traps)))
trap.list <- melt(trap.list)
names(trap.list)<-c("traps","site")
full.com<-left_join(full.com,trap.list,join_by("site"=="site"))

### Plot observed vs Estimated with errors around variance

# Determine overlap vs. no

full.com$signif.rich<- ifelse(full.com$Observed.rich > (full.com$Estimator.rich + full.com$Est_s.e..rich) |
                           full.com$Observed.rich < (full.com$Estimator.rich - full.com$Est_s.e..rich), 
                         1, 0)

full.com$signif.div<- ifelse(full.com$Observed.div > (full.com$Estimator.div + full.com$Est_s.e..div) |
                                full.com$Observed.div < (full.com$Estimator.div - full.com$Est_s.e..div), 
                              1, 0)

# Create rank order observed richness values
full.com<-full.com[order(full.com$Observed.rich,decreasing=TRUE),]
full.com$obsRank.rich<-1:nrow(full.com)

# Model proportion of estimated hill number observed to date as a function of the estimated value, turnover and number of years of sampling

summary(glm(propObs.rich~Estimator.rich+turnover+years,family='quasibinomial',full.com))
summary(glm(propObs.div~Estimator.div+turnover+years,family='quasibinomial',full.com))

# Plot the richness and turnover relationships

rich <- ggplot(full.com, aes(x = turnover, y = propObs.rich, color = as.numeric(Estimator.rich))) +
  geom_point(aes(size = years)) +
  geom_smooth(method = "glm", method.args = list(family = "quasibinomial"), color = "black") +
  labs(x = "", y = "Observed/Estimated Richness") +
  theme_minimal() +
  annotate("text", x = -Inf, y = Inf, label = "a", fontface = "bold", hjust = -0.2, vjust = 1.3, size = 6) +
  scale_color_viridis_c(option = "D", name = "Est. richness")

div <- ggplot(full.com, aes(x = turnover, y = propObs.div, color = as.numeric(Estimator.div))) +
  geom_point(aes(size = years)) +
  geom_smooth(method = "glm", method.args = list(family = "quasibinomial"), color = "black") +
  labs(x = "Mean Species Turnover", y = "Observed/Estimated Diversity") +
  theme_minimal() +
  annotate("text", x = -Inf, y = Inf, label = "b", fontface = "bold", hjust = -0.1, vjust = 3, size = 6) +
  scale_color_viridis_c(option = "D", name = "Est. diversity")

combined <- rich / div + 
  theme(legend.position = "right")

print(combined)

# Plot overlap in observed vs estimated Hill numbers

rich.plot <- 
  ggplot(full.com, 
         aes(x = obsRank.rich, y = Observed.rich)) +
  geom_rect(aes(xmin = obsRank.rich - 0.5, 
                xmax = obsRank.rich + 0.5, 
                ymin = -Inf, 
                ymax = Inf, alpha = signif.rich),
            fill = "grey") +
  geom_errorbar(aes(ymin = Estimator.rich - Est_s.e..rich, 
                    ymax = Estimator.rich + Est_s.e..rich, 
                    color = turnover), width = 0, size = 2) + 
  geom_point(color = "black", size = 2) + 
  labs(x = "", y = "Richness") +
  annotate("text", x = -Inf, y = Inf, label = "a", fontface = "bold", hjust = -0.2, vjust = 1.3, size = 6) +
  theme_bw() + 
  scale_color_viridis_c(option = "D", name = "turnover") +
  scale_alpha(range = c(0, 0.5), guide = "none") +
  theme(
    axis.title.y = element_text(size = 16),  
    axis.text.y = element_text(size = 14),   
    axis.title.x = element_blank(),  
    axis.text.x = element_blank(),  
    legend.position = "none",
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    panel.spacing = unit(0, "cm")
  )

div.plot <- 
  ggplot(full.com, 
         aes(x = obsRank.rich, y = Observed.div)) +
  geom_rect(aes(xmin = obsRank.rich - 0.5, 
                xmax = obsRank.rich + 0.5, 
                ymin = -Inf, 
                ymax = Inf, alpha = signif.div),
            fill = "grey") +
  geom_errorbar(aes(ymin = Estimator.div - Est_s.e..div, 
                    ymax = Estimator.div + Est_s.e..div, 
                    color = turnover), width = 0, size = 2) + 
  geom_point(color = "black", size = 2) + 
  labs(x = "Rank-order Observed Richness", y = "Diversity") +
  annotate("text", x = -Inf, y = Inf, label = "b", fontface = "bold", hjust = -0.2, vjust = 1.3, size = 6) +
  theme_bw() + 
  scale_color_viridis_c(option = "D", name = "turnover") +
  scale_alpha(range = c(0, 0.5), guide = "none") +
  theme(
    axis.title = element_text(size = 16),  
    axis.text = element_text(size = 14),
    legend.position = c(1, 1),  # Adjusted to move inside the plot
    legend.justification = c(1, 1), # Anchor the legend to the top right
    legend.background = element_rect(fill = NA, color = NA),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    panel.spacing = unit(0, "cm")
  )

plot_grid(rich.plot, div.plot, ncol = 1, align = 'v', axis = 'tb')


# Violin plots of turnover vs. overlap with estimator +/- se
full.com$signifText.rich<-"Overlap"
full.com$signifText.div<-"Overlap"
full.com$signifText.rich[which(full.com$signif.rich==1)]<-"No Overlap"
full.com$signifText.div[which(full.com$signif.div==1)]<-"No Overlap"

richplot<-ggplot(full.com,aes(x=turnover, y=as.factor(signifText.rich)))+
  labs(y = "Richness", x = "Mean Species Turnover") +
  geom_violin() + stat_summary(
    fun = "mean", geom = "point", shape = 20, size = 3, color = "black")+
  theme_minimal()   + theme(
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(),    
    axis.title.y = element_text(size = 18),  
    axis.text.y = element_text(size = 16)     
  )

divplot<-ggplot(full.com,aes(x=turnover, y=as.factor(signifText.div)))+
  labs(y = "Diversity", x = "Mean Species Turnover") +
  geom_violin() + stat_summary(
    fun = "mean", geom = "point", shape = 20, size = 3, color = "black")+
  theme_minimal()   + theme(
    axis.title = element_text(size = 18),  
    axis.text = element_text(size = 16)  
  )

combined_plot <- richplot/divplot
combined_plot


signif_model_rich<-glm(signif.rich~turnover,full.com,family="binomial")
summary(signif_model_rich)

signif_model_div<-glm(signif.div~turnover,full.com,family="binomial")
summary(signif_model_div)


