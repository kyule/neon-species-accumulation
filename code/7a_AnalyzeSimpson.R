### This is the primary file for formatting and plotting NEON species accumulation data

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
library("patchwork")
library("cowplot")

# Load in the formatted clean data
#Make sure the results are correctly configured

load(file=paste0(datapath,"iNEXTandTurnoverResults_withSimpson.Robj"))

### Pull Estimated Asymptotic and Observed simpson and Diversity Values out of the results list -- full data only

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

# Extract Simpson
Simpson <- asyest[grep("Simpson", rownames(asyest)), ]
Simpson <- Simpson[, c(1:4, ncol(Simpson))]  

full.com<-Simpson

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
full.com$propObs<-full.com$Observed/full.com$Estimator

# Remove GUAN due to very low data availability
full.com<-full.com[which(full.com$site!="GUAN"),]

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

full.com$final.est<-full.com$prop.final.est<-full.com$final.est.se<-full.com$final.obs<-full.com$final.obs.div<-0

for (i in 1:nrow(full.com)){
  sim.est<-full.com$Estimator[which(full.com$site==full.com$site[i] & full.com$year=="full")]

  full.com$final.est.se[i]<-full.com$Est_s.e.[which(full.com$site==full.com$site[i] & full.com$year=="full")]

  full.com$final.est[i]<-sim.est

  full.com$prop.final.est[i]<-full.com$Observed[i]/sim.est

  full.com$final.obs[i]<-full.com$Observed[which(full.com$site==full.com$site[i] & full.com$year=="full")]
  
}

## Add warning or not
names(results[["warnings"]])[2]<-'warning'
full.com<-left_join(full.com,results[["warnings"]],join_by("site"=="site","year"=="year"))
full.com$warning[which(is.na(full.com$warning)==FALSE)]<-"Y"
full.com$warning[which(is.na(full.com$warning)==TRUE)]<-"N"

# Determine overlap vs. no

full.com$signif<- ifelse(full.com$Observed > (full.com$Estimator + full.com$Est_s.e.) |
                                full.com$Observed < (full.com$Estimator - full.com$Est_s.e.), 
                              1, 0)

full.com$final.signif<- ifelse(full.com$Observed > (full.com$final.est + full.com$final.est.se) |
                                full.com$Observed < (full.com$final.est - full.com$final.est.se), 
                              1, 0)


full.com$signifText<-"Overlap"
full.com$signifText[which(full.com$signif==1)]<-"No Overlap"

# Create rank order observed simpson values
full<-full.com[which(full.com$year=="full"),]
full<-full[order(full$Observed,decreasing=TRUE),]
full$obsRank<-1:nrow(full)

full.com$obsRank<-0

for (i in 1:nrow(full)){
  full.com$obsRank[which(full.com$site==full$site[i])]<-full$obsRank[i]
}


full<- full.com[which(full.com$year=="full"),]




# Plot the simpson and turnover relationships

sim_limits <- range(full.com$final.est, na.rm = TRUE)

sim.all <- ggplot(full.com[which((full.com$year != "full") & !is.na(full.com$turnover)),], aes(x = turnover, y = prop.final.est, color = as.numeric(Estimator))) +
  geom_point(aes(size = years)) +
  scale_size_continuous(range = c(1, 3)) +
  geom_smooth(method = "glm", color = "black") +
  labs(x = "Years of sampling", y = "(Observed Simpson Diversity)/(Estimated Simpson Diversity)") +
  theme_minimal() +
  #annotate("text", x = -Inf, y = Inf, label = "a", fontface = "bold", hjust = -0.2, vjust = 1.3, size = 6) +
  scale_color_viridis_c(option = "D", name = "Est. Simpson Diversity",limits = sim_limits)

print(sim.all)


# Set shared color limits
turnover_limits <- range(c(full$turnover, full.com$turnover), na.rm = TRUE)

# simpson
sim.plot <- 
  ggplot(full, aes(x = obsRank, y = Observed)) +
  geom_errorbar(aes(ymin = Estimator - Est_s.e., 
                    ymax = Estimator + Est_s.e., 
                    color = turnover), width = 0, size = 2) + 
  geom_point(color = "black", size = 2) + 
  labs(x = "Rank order Simpson Diversity", y = "Simpson Diversity") +
  theme_bw() + 
  scale_color_viridis_c(option = "D", name = "Mean Turnover", limits = turnover_limits) +
  scale_alpha(range = c(0, 0.5), guide = "none") +
  theme(
    axis.title = element_text(size = 12),  
    axis.text = element_text(size = 12),
    legend.position = c(1, 1),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = NA, color = NA),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    panel.spacing = unit(0, "cm")
  )

print(sim.plot)

# simpson Proportion
sim.prop <- ggplot(full.com, aes(x = years, y = prop.final.est, color = turnover)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(x = "Years of sampling", y = "(Observed Simpson Diversity)/(Estimated Simpson Diversity)") +
  geom_smooth(color = 'black') +
  theme(legend.position = "none") +
  theme(
    axis.title = element_text(size = 12),  
    axis.text = element_text(size = 12),
    legend.position = "none",
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    panel.spacing = unit(0, "cm")
  ) +
  scale_color_viridis_c(option = "D", name = "Turnover", limits = turnover_limits)
print(sim.prop)

#write.csv(full.com,paste0(datapath,'communityResults_withSimpson.csv'))
