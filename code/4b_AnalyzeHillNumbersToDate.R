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

load(file=paste0(datapath,"iNEXTandTurnoverResults.Robj"))

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

full.com$final.est.rich<-full.com$final.est.div<-full.com$prop.final.est.rich<-full.com$prop.final.est.div<-full.com$final.est.rich.se<-full.com$final.est.div.se<-full.com$final.obs.rich<-full.com$final.obs.div<-0

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
  full.com$final.obs.div[i]<-full.com$Observed.div[which(full.com$site==full.com$site[i] & full.com$year=="full")]
  
  
}

## Add warning or not
names(results[["warnings"]])[2]<-'warning'
full.com<-left_join(full.com,results[["warnings"]],join_by("site"=="site","year"=="year"))
full.com$warning[which(is.na(full.com$warning)==FALSE)]<-"Y"
full.com$warning[which(is.na(full.com$warning)==TRUE)]<-"N"

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

full.com$signifText.rich<-"Overlap"
full.com$signifText.div<-"Overlap"
full.com$signifText.rich[which(full.com$signif.rich==1)]<-"No Overlap"
full.com$signifText.div[which(full.com$signif.div==1)]<-"No Overlap"

# Create rank order observed richness values
full<-full.com[which(full.com$year=="full"),]
full<-full[order(full$Observed.rich,decreasing=TRUE),]
full$obsRank<-1:nrow(full)

full.com$obsRank<-0

for (i in 1:nrow(full)){
  full.com$obsRank[which(full.com$site==full$site[i])]<-full$obsRank[i]
}

# Remove GUAN due to very low data availability
full.com<-full.com[which(full.com$site!="GUAN"),]
full<- full.com[which(full.com$year=="full"),]




# Plot the richness and turnover relationships

rich_limits <- range(full.com$final.est.rich, na.rm = TRUE)
div_limits <- range(full.com$final.est.div, na.rm = TRUE)

rich.all <- ggplot(
  full.com[which((full.com$year != "full") & !is.na(full.com$turnover)),],
  aes(x = turnover, y = prop.final.est.rich, color = as.numeric(Estimator.rich))
) +
  scale_color_viridis_c(
    option = "D",
    name = "Estimated richness",
    limits = rich_limits,
    guide = guide_colorbar(order = 1)   # <-- puts this legend first (on top)
  ) +
  geom_point(aes(size = years)) +
  scale_size_continuous(
    range = c(1, 3),
    guide = guide_legend(order = 2)      # <-- puts size legend below
  ) +
  geom_smooth(method = "glm", method.args = list(family = "quasibinomial"), color = "black") +
  labs(x = "", y = "(Observed richness)/(estimated richness)") +
  theme_minimal() +
  annotate("text", x = 0.22, y = 1, label = "a)", fontface = "bold")


div.all <- ggplot(
  full.com[which((full.com$year != "full") & !is.na(full.com$turnover)),],
  aes(x = turnover, y = prop.final.est.div, color = as.numeric(Estimator.div))
) +
  geom_point(aes(size = years)) +
  scale_size_continuous(
    range = c(1, 3),
    guide = guide_legend(order = 2)   # only matters if you later show it
  ) +
  geom_smooth(method = "lm", color = "black") +
  labs(x = "Mean species turnover", y = "(Observed diversity)/(estimated diversity)") +
  theme_minimal() +
  annotate("text", x = 0.22, y = 1.15, label = "b)", fontface = "bold") +
  scale_color_viridis_c(
    option = "D",
    name = "Estimated diversity",
    limits = div_limits,
    guide = guide_colorbar(order = 1)  # <-- puts color scale bar first
  ) 


combined.all <- rich.all / div.all + 
  theme(legend.position = "right")

print(combined.all)


# Set shared color limits
turnover_limits <- range(c(full$turnover, full.com$turnover), na.rm = TRUE)

# Plot 1: Richness
rich.plot <- 
  ggplot(full, aes(x = obsRank, y = Observed.rich)) +
  geom_rect(aes(xmin = obsRank - 0.5, xmax = obsRank + 0.5, ymin = -Inf, ymax = Inf, alpha = signif.rich),
            fill = "grey") +
  geom_errorbar(aes(ymin = Estimator.rich - Est_s.e..rich, 
                    ymax = Estimator.rich + Est_s.e..rich, 
                    color = turnover), width = 0, size = 2) + 
  geom_point(color = "black", size = 2) + 
  labs(x = "", y = "Richness") +
  annotate("text", x = 0, y = 300, label = "a)", size = 5, fontface = "bold") +
  theme_bw() + 
  scale_color_viridis_c(option = "D", name = "Mean turnover", limits = turnover_limits) +
  scale_alpha(range = c(0, 0.5), guide = "none") +
  theme(
    axis.title = element_text(size = 10),  
    axis.text = element_text(size = 10),
    legend.position = c(1, 1),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = NA, color = NA),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    subtitle ='a',
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    panel.spacing = unit(0, "cm")
  )


# Plot 2: Diversity
div.plot <- 
  ggplot(full, aes(x = obsRank, y = Observed.div)) +
  geom_rect(aes(xmin = obsRank - 0.5, xmax = obsRank + 0.5, ymin = -Inf, ymax = Inf, alpha = signif.div),
            fill = "grey") +
  geom_errorbar(aes(ymin = Estimator.div - Est_s.e..div, 
                    ymax = Estimator.div + Est_s.e..div, 
                    color = turnover), width = 0, size = 2) + 
  geom_point(color = "black", size = 2) + 
  labs(x = "Rank-order observed Richness", y = "Shannon diversity") +
  annotate("text", x = 0, y = 40, label = "b)", size = 5, fontface = "bold") +
  theme_bw() + 
  scale_color_viridis_c(option = "D", name = "Turnover", limits = turnover_limits) +
  scale_alpha(range = c(0, 0.5), guide = "none") +
  theme(
    axis.title.y = element_text(size = 10),  
    axis.text.y = element_text(size = 10),   
    axis.title.x = element_text(size = 10),  
    axis.text.x = element_text(size = 10),  
    legend.position = "none",
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    panel.spacing = unit(0, "cm")
  )

# Combine first pair
hill_combo <- plot_grid(rich.plot, div.plot, ncol = 1, align = 'v', axis = 'tb')

# Plot 3: Richness Proportion
rich.prop <- ggplot(full.com, aes(x = years, y = prop.final.est.rich, color = turnover)) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = c(2, 4,6, 8))+
  theme_minimal() +
  labs(x = "", y = "(Observed richness)/(estimated richness)") +
  geom_smooth(color = 'black') +
  annotate("text", x = 1, y = 1, label = "c)", size = 5, fontface = "bold") +
  theme(legend.position = "none") +
  theme(
    axis.title = element_text(size = 10),  
    axis.text = element_text(size = 10),
    legend.position = "none",
    plot.margin = unit(c(0.15, 0, 0, 0), "cm"),
    panel.spacing = unit(0, "cm")
  ) +
  scale_color_viridis_c(option = "D", name = "Turnover", limits = turnover_limits)

# Plot 4: Diversity Proportion
div.prop <- ggplot(full.com, aes(x = years, y = prop.final.est.div, color = turnover)) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = c(2, 4,6, 8))+
  theme_minimal() +
  labs(x = "Years of sampling", y = "(Observed diversity)/(estimated diversity)") +
  annotate("text", x = 1, y = 1.35, label = "d)", size = 5, fontface = "bold") +
  geom_smooth(color = 'black') +
  theme(legend.position = "none") +
  theme(
    axis.title.y = element_text(size = 10),  
    axis.text.y = element_text(size = 10),   
    axis.title.x = element_text(size = 10),  
    axis.text.x = element_text(size = 10),  
    legend.position = "none",
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    panel.spacing = unit(0, "cm")
  ) +
  scale_color_viridis_c(option = "D", name = "Turnover", limits = turnover_limits)

prop_combo <- plot_grid(rich.prop, div.prop, ncol = 1, align = 'v', axis = 'lr', rel_heights = c(1, 1))


final_combo <- plot_grid(hill_combo, plot_spacer()+ theme_void(), prop_combo, 
                         ncol = 3, 
                         rel_widths = c(1, 0.05, 0.6), 
                         align = 'h')
tiff(paste0(datapath,"figure1.tiff"), units="px", width=1961, height=1700, res=300)
final_combo
dev.off()

