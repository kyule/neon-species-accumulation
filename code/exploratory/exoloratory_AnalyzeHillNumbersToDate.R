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
library("purrr")
library("glmtoolbox")
library('lme4')

# Load in the formatted clean data, or download and create it. 
#Make sure the results are correctly configured

load(file=paste0(datapath,"iNEXTandTurnoverResults.Robj"))

# Remove GUAN due to very low number of beetles captured. 
### Only 36 beetles in 7 years of sampling, almost complete turnover in communities between most years
#results<- results[!names(results)=="GUAN"]

### Pull Estimated Asymptotic and Observed Richness and Diversity Values out of the results list -- full data only

asyest.list <- lapply(results, function(site) {
  bind_rows(
    lapply(names(site), function(year) {  # Iterate through all result years
      if (is.list(site[[year]]) && !is.null(site[[year]]$out$AsyEst)) {  # Check if AsyEst exists
        df <- data.frame(site[[year]]$out$AsyEst)
        df$year <- year
        return(df)
      } else {
        return(NULL)  # Skip missing values
      }
    }),
    .id = "site"
  )
})

# Combine results into one data frame
asyest <- bind_rows(asyest.list, .id = "site")

# Extract Shannon Diversity
Diversity <- asyest[grep("Shannon", rownames(asyest)), ]
Diversity <- Diversity[, c(1:4, ncol(Diversity))]  

# Extract Richness
Richness <- asyest[grep("Richness", rownames(asyest)), ]
Richness <- Richness[, c(1:4, ncol(Richness))]  


full.com <- left_join(Richness,Diversity,join_by("site"=="site","year"=="year"),suffix = c(".rich",".div"))

#### Pull the turnovers out of the results list and input into the full community data frame
# Calculate cumulative mean for each site and year
turnover.list <- lapply(results, function(item) data.frame(item$turnover))
turnover.list <- bind_rows(turnover.list, .id = "site")

turnover <- turnover.list %>%
  group_by(site) %>%
  arrange(year) %>%
  mutate(turnover = cummean(total))  # This gives the cumulative mean for each site

# Calculate the mean for all years for each site and add a row for 'full'
turnover_full <- turnover.list %>%  # Use the original data, not the cumulative mean
  group_by(site) %>%
  summarise(year = "full", turnover = mean(total))  # Mean of total, not turnover

# Ensure year is a character for consistency
turnover$year <- as.character(turnover$year)

# Bind the full dataset
turnover <- bind_rows(turnover, turnover_full)
turnover<-turnover[,which(names(turnover) %in% c("site","year","turnover"))]

# Bind turnover to the full community dataset

full.com <- left_join(full.com,turnover,join_by("site"=="site","year"=="year"))

full.com$propObs.rich<-full.com$Observed.rich/full.com$Estimator.rich
full.com$propObs.div<-full.com$Observed.div/full.com$Estimator.div


# Take into account the number of traps that have been sampled
trap.df <- bind_rows(lapply(names(results), function(site_name) {  # Iterate over sites
  site <- results[[site_name]]
  
  bind_rows(lapply(names(site), function(year) {  # Iterate over years
    if (is.list(site[[year]]) && !is.null(site[[year]]$traps)) {  
      data.frame(
        site = site_name,
        year = year,
        traps = length(site[[year]]$traps)  # Get the length of the traps vector
      )
    } else {
      return(NULL)  # Skip missing values
    }
  }))
}))

full.com<-left_join(full.com,trap.df,join_by("site"=="site","year"=="year"))

# Add number of years of sampling

full.com <- full.com %>%
  group_by(site) %>%             
  arrange(year, .by_group = TRUE) %>% 
  mutate(years = row_number())

full.com$years[which(full.com$year=="full")]<-full.com$years[which(full.com$year=="full")]-1

# Define most complete estimators and associated proportions

full.com$final.est.rich<-full.com$final.est.div<-full.com$prop.final.est.rich<-full.com$prop.final.est.div<-full.com$final.est.rich.se<-full.com$final.est.div.se<-full.com$final.obs.rich<-0

for (i in 1:nrow(full.com)){
  rich.est<-full.com$Estimator.rich[which(full.com$site==full.com$site[i] & full.com$year=="full")]
  div.est<-full.com$Estimator.div[which(full.com$site==full.com$site[i] & full.com$year=="full")]
  
  full.com$final.est.rich.se[i]<-full.com$Est_s.e..rich[which(full.com$site==full.com$site[i] & full.com$year=="full")]
  full.com$final.est.div.se[i]<-full.com$Est_s.e..div[which(full.com$site==full.com$site[i] & full.com$year=="full")]
  
  full.com$final.est.rich[i]<-rich.est
  full.com$final.est.div[i]<-div.est
  
  full.com$prop.final.est.rich[i]<-full.com$Observed.rich[i]/rich.est
  full.com$prop.final.est.div[i]<-full.com$Observed.div[i]/div.est
  
  full.com$final.obs.rich[i]<-full.com$Observed.rich[which(full.com$site==full.com$site[i] & full.com$year=="full")]
  
  
}

## Add warning or not
names(results[["warnings"]])[2]<-'warning'
full.com<-left_join(full.com,results[["warnings"]],join_by("site"=="site","year"=="year"))
full.com$warning[which(is.na(full.com$warning)==FALSE)]<-"Y"
full.com$warning[which(is.na(full.com$warning)==TRUE)]<-"N"


### Plot observed vs Estimated with errors around variance

# Determine overlap vs. no

full.com$signif.rich<- ifelse(full.com$Observed.rich > (full.com$Estimator.rich + full.com$Est_s.e..rich) |
                                full.com$Observed.rich < (full.com$Estimator.rich - full.com$Est_s.e..rich), 
                              1, 0)

full.com$signif.div<- ifelse(full.com$Observed.div > (full.com$Estimator.div + full.com$Est_s.e..div) |
                               full.com$Observed.div < (full.com$Estimator.div - full.com$Est_s.e..div), 
                             1, 0)

full.com$final.signif.rich<- ifelse(full.com$Observed.rich > (full.com$final.est.rich + full.com$final.est.rich.se) |
                                full.com$Observed.rich < (full.com$final.est.rich - full.com$final.est.rich.se), 
                              1, 0)

full.com$final.signif.div<- ifelse(full.com$Observed.div > (full.com$final.est.div + full.com$final.est.div.se) |
                               full.com$Observed.div < (full.com$final.est.div + full.com$final.est.div.se), 
                             1, 0)

# Create rank order observed richness values
full<-full.com[which(full.com$year=="full"),]
full<-full[order(full$Observed.rich,decreasing=TRUE),]
full$obsRank<-1:nrow(full)

full.com$obsRank<-0

for (i in 1:nrow(full)){
  full.com$obsRank[which(full.com$site==full$site[i])]<-full$obsRank[i]
}



###### Model proportion of estimated hill number observed to date as a function of the estimated value, turnover and number of years of sampling

a<-glm(propObs.rich~ turnover*Estimator.rich*years,family='quasibinomial',full)
summary(a)
stepCriterion(a)
b<-glm(propObs.rich~ Estimator.rich*years,family='quasibinomial',full)
summary(b)

a<-glm(propObs.div~turnover*Estimator.div*years,family='quasibinomial',full)
summary(a)
stepCriterion(a)
b<-glm(propObs.div~turnover,family='quasibinomial',full)
summary(b)

# Plot proportion by number of years

ggplot(full.com,aes(x=years,y=prop.final.est.rich))+
          geom_point() +
         geom_smooth()
ggplot(full.com,aes(x=years,y=prop.final.est.div))+
  geom_point() +
  geom_smooth()

# Plot the richness and turnover relationships

rich <- ggplot(full, aes(x = turnover, y = prop.final.est.rich, color = as.numeric(Estimator.rich))) +
  geom_point(aes(size = years)) +
  geom_smooth(method = "glm", method.args = list(family = "quasibinomial"), color = "black") +
  labs(x = "", y = "Observed/Estimated Richness") +
  theme_minimal() +
  annotate("text", x = -Inf, y = Inf, label = "a", fontface = "bold", hjust = -0.2, vjust = 1.3, size = 6) +
  scale_color_viridis_c(option = "D", name = "Est. richness")

div <- ggplot(full, aes(x = turnover, y = prop.final.est.div, color = as.numeric(Estimator.div))) +
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
  ggplot(full.com, 
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


