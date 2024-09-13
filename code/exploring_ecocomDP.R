## ----load libraries, eval=F, comment=NA-----------------------------------------------------------------------------------------------

# clean out workspace

rm(list = ls()) # OPTIONAL - clear out your environment
gc()            # Uncomment these lines if desired

# load packages
library(tidyverse)
library(neonUtilities)
library(vegan)
library(ecocomDP)
library(Hmisc)

# source .r file with my NEON_TOKEN
source("/Users/kelsey/Github/neon-biorepo-tools/configini.R")

## ----download-macroinvert, message=FALSE, warning=FALSE, results='hide'---------------------------------------------------------------

# search for invertebrate data products
my_search_result <- 
    ecocomDP::search_data(text = "invertebrate")
View(my_search_result)

# pull data for the NEON  "beetle
# collection"
my_data <- ecocomDP::read_data(
  id = "neon.ecocomdp.10022.001.001",
  site = "all",
  startdate = "2020-01",
  enddate = "2024-07",
  token = Neon_Token,
  check.size = FALSE)



## ----view-ecocomDP-str----------------------------------------------------------------------------------------------------------------

# examine the structure of the data object that is returned
my_data %>% names()
my_data$id

# short list of package summary data
my_data$metadata$data_package_info

# validation issues? None if returns an empty list
my_data$validation_issues

# examine the tables
my_data$tables %>% names()
my_data$tables$taxon %>% head()
my_data$tables$observation %>% head()


# Explore the spatial and temporal coverage 
# of the dataset
my_data %>% ecocomDP::plot_sample_space_time()


# Explore the taxonomic resolution in the dataset. 
# What is the most common taxonomic resolution (rank) 
# for macroinvertebrate identifications in this dataset?
my_data %>% ecocomDP::plot_taxa_rank()

## ----flattening-and-cleaning, message=FALSE, warning=FALSE----------------------------------------------------------------------------

# flatten the ecocomDP data tables into one flat table
my_data_flat <- my_data$tables %>% ecocomDP::flatten_data()



## ----data-vis-richness-time, message=FALSE, warning=FALSE, fig.cap="Benthic algal richness by year at ARIK and COMO"------------

# plot richness by year for a site
my_data_flat %>% filter(siteID =="SRER") %>% ecocomDP::plot_taxa_diversity(time_window_size = "year")

## ----more-cleaning, message=FALSE, warning=FALSE--------------------------------------------------------------------------------------

# Note that for this data product
# neon_sample_id = event_id
# event_id is a grouping variable for the observation 
# table in the ecocomDP data model

# Check for multiple taxon counts per taxon_id by 
# event_id. 
summary_data<-my_data_flat %>% 
                  mutate(year=year(datetime)) %>% 
                  mutate(site_year=paste(siteID,year,sep=".")) %>% 
                  group_by(site_year,taxon_id) %>%
                  summarize(abundance = sum(value*as.numeric(trappingDays))) %>%
                  dplyr::filter(abundance > 1)


## ----Species accumulation plot. Confidence intervals are based on random permutations of observed samples."----

# make data wide for a single site
my_data_wide <- summary_data %>% 
  filter(str_detect(site_year, "ABBY")) %>%
  pivot_wider(names_from = taxon_id, 
              values_from = abundance,
              values_fill =  0) %>%
  tibble::column_to_rownames("site_year")
  
# Calculate and plot species accumulcation curve for a site
# The CIs are based on random permutations of observed samples
spec_accum_result <- my_data_wide %>% 
                          vegan::specaccum(., "random")
plot(spec_accum_result)



## ----compare-obs-sim-SAC--------------------------------------------------------------------------------------------------------------

# Extract the resampling data used in the above algorithm
spec_resamp_data <- data.frame(
  data_set = "observed", 
  sampling_effort = rep(1:nrow(spec_accum_result$perm), 
                        each = ncol(spec_accum_result$perm)),
  richness = c(t(spec_accum_result$perm)))


# Fit species accumulation model
spec_accum_mod_1 <- my_data_wide %>% vegan::fitspecaccum(model = "arrh")


# create a "predicted" data set from the model to extrapolate out 
# beyond the number of samples collected
sim_spec_data <- data.frame()
for(i in 1:27){
  d_tmp <- data.frame(
    data_set = "predicted",
    sampling_effort = i,
    richness = predict(spec_accum_mod_1, newdata = i))
  
  sim_spec_data <- sim_spec_data %>%
    bind_rows(d_tmp)
}


# plot the "observed" and "simulated" curves with 95% CIs
data_plot <- spec_resamp_data %>% bind_rows(sim_spec_data) 

data_plot %>%
  ggplot(aes(sampling_effort, richness, 
             color = as.factor(data_set),
             fill = as.factor(data_set),
             linetype = as.factor(data_set))) +
  stat_summary(fun.data = median_hilow, fun.args = list(conf.int = .95), 
               geom = "ribbon", alpha = 0.25) +
  stat_summary(fun.data = median_hilow, geom = "line", 
               size = 1) 
    

