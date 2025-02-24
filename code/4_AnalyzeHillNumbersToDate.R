### This is the primary file for plotting NEON species accumulation data

#### BEFORE RUNNING
# Users must define their own paths

# datapath<-"user defined path"
# codepath<-"user defined path"

# And users should remove below line that loads in personal paths
source("/Users/kelsey/Github/neon-species-accumulation/configini.R")

#set seed
set.seed(85705)

# Load some necessary packages

library("ggplot2")
library("patchwork")
library("cowplot")

# Load in the formatted clean data
#Make sure the results are correctly configured

load(file=paste0(datapath,"iNEXTandTurnoverResults.Robj"))

# Plot proportion by number of years

rich.prop<-ggplot(full.com,aes(x=years,y=prop.final.est.rich,color=turnover))+
          geom_point(size=3) +
          theme_minimal() +
          labs(x = "", y = "Observed/Estimated Richness") +
          geom_smooth(color='black') +
          scale_color_viridis_c(option = "D", name = "turnover")

div.prop<-ggplot(full.com,aes(x=years,y=prop.final.est.div,color=turnover))+
  geom_point(size=3) +
  theme_minimal() +
  labs(x = "Years of sampling", y = "Observed/Estimated Diversity") +
  geom_smooth(color='black')+
  scale_color_viridis_c(option = "D", name = "turnover")

combo <- rich.prop / div.prop + 
  theme(legend.position = "right")

print(combo)


# Plot the richness and turnover relationships

rich <- ggplot(full.com[which((full.com$year != "full") & !is.na(full.com$turnover)),], aes(x = turnover, y = prop.final.est.rich, color = as.numeric(Estimator.rich))) +
  geom_point(aes(size = years)) +
  geom_smooth(method = "glm", method.args = list(family = "quasibinomial"), color = "black") +
  labs(x = "", y = "Observed/Estimated Richness") +
  theme_minimal() +
  annotate("text", x = -Inf, y = Inf, label = "a", fontface = "bold", hjust = -0.2, vjust = 1.3, size = 6) +
  scale_color_viridis_c(option = "D", name = "Est. richness")

div <- ggplot(full.com[which((full.com$year != "full") & !is.na(full.com$turnover)),], aes(x = turnover, y = prop.final.est.div, color = as.numeric(Estimator.div))) +
  geom_point(aes(size = years)) +
  geom_smooth(method = "lm", color = "black") +
  labs(x = "Mean Species Turnover", y = "Observed/Estimated Diversity") +
  theme_minimal() +
  annotate("text", x = -Inf, y = Inf, label = "b", fontface = "bold", hjust = -0.1, vjust = 3, size = 6) +
  scale_color_viridis_c(option = "D", name = "Est. diversity")

combined <- rich / div + 
  theme(legend.position = "right")

print(combined)

# Plot overlap in observed vs estimated Hill numbers

rich.plot <- 
  ggplot(full, 
         aes(x = obsRank, y = Observed.rich)) +
  geom_rect(aes(xmin = obsRank - 0.5, 
                xmax = obsRank + 0.5, 
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
  ggplot(full, 
         aes(x = obsRank, y = Observed.div)) +
  geom_rect(aes(xmin = obsRank - 0.5, 
                xmax = obsRank + 0.5, 
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


richplot<-ggplot(full,aes(x=turnover, y=as.factor(signifText.rich)))+
  labs(y = "Richness", x = "Mean Species Turnover") +
  geom_violin() + stat_summary(
    fun = "mean", geom = "point", shape = 20, size = 3, color = "black")+
  theme_minimal()   + theme(
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(),    
    axis.title.y = element_text(size = 18),  
    axis.text.y = element_text(size = 16)     
  )

divplot<-ggplot(full,aes(x=turnover, y=as.factor(signifText.div)))+
  labs(y = "Diversity", x = "Mean Species Turnover") +
  geom_violin() + stat_summary(
    fun = "mean", geom = "point", shape = 20, size = 3, color = "black")+
  theme_minimal()   + theme(
    axis.title = element_text(size = 18),  
    axis.text = element_text(size = 16)  
  )

combined_plot <- richplot/divplot
combined_plot


### save the results file for stat analysis

write.csv(full.com,paste0(datapath,'communityResults.csv'),row.names=FALSE)




