### This is the primary file for formatting NEON species accumulation data

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


