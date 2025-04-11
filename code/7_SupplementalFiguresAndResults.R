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
library('dplyr')
library('lme4')
library("glmtoolbox")
library('lmerTest')
library('gt')
library('gtsummary')

# Load in the formatted clean data
# Make sure the results are correctly configured

full.com<-read.csv(paste0(datapath,'communityResults.csv'))
rich.thresh<-read.csv(paste0(datapath,'richnessThresh.csv'))
div.thresh<-read.csv(paste0(datapath,'diversityThresh.csv'))

# Make richness and diversity table

fullcom.summary <- full.com %>%
  filter(year == "full") %>%
  select(site, final.obs.rich, final.est.rich, final.est.rich.se, 
         final.obs.div, final.est.div, final.est.div.se, turnover, years) %>%
  mutate(
    final.est.rich_combined = sprintf("%.2f (±%.2f)", final.est.rich, final.est.rich.se),
    final.est.div_combined = sprintf("%.2f (±%.2f)", final.est.div, final.est.div.se)
  ) %>%
  select(-final.est.rich, -final.est.rich.se, -final.est.div, -final.est.div.se)  

gt_table<-gt(fullcom.summary) %>%
  fmt_number(columns = 'final.obs.rich', decimals = 0) %>%
  fmt_number(columns = "final.obs.div", decimals = 2,drop_trailing_zeros = FALSE) %>%
  fmt_number(columns = "turnover", decimals = 2,drop_trailing_zeros = FALSE) %>%
  fmt_number(columns = "years", decimals = 0) %>%
  cols_label(
    site = "Site",
    final.obs.rich = "Observed",
    final.est.rich_combined = "Estimated (±SE)",
    final.obs.div = "Observed",
    final.est.div_combined = "Estimated (±SE)",
    turnover = "Mean Turnover",
    years = "Sampling Years"
  ) %>%
  tab_spanner(
    label = "Richness",
    columns = c(final.obs.rich, final.est.rich_combined)
  ) %>%
  tab_spanner(
    label = "Diversity",
    columns = c(final.obs.div, final.est.div_combined)
  ) %>%
  cols_align(
    align = "right",
    columns = c(final.est.rich_combined, final.est.div_combined) 
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything()) 
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_spanners(everything())  
  ) %>%
  tab_style(
    style = cell_fill(color = "#f5f5f5"), 
    locations = cells_body(rows = seq(2, nrow(fullcom.summary), by = 2))  
  )

gtsave(gt_table, filename = "fullcom_table.png", path = datapath)

### Turnover by overlap

signif_model_rich<-glmer(signif.rich~ turnover
                         + (1|site),
                         data=mod_dat,
                         family="binomial")
summary(signif_model_rich)

# prob of overlap decreases with turnover

signif_model_div<-glmer(signif.div~ turnover
                        + (1|site),
                        data=mod_dat,
                        family="binomial")
summary(signif_model_div)

# prob of overlap decreases with turnover


#### Full data only analysis

mod_dat<-full.com[which(full.com$year=="full"),] # full data only

mod_dat$years<-scale(mod_dat$years,center=TRUE,scale=TRUE)
mod_dat$fin.est.rich<-scale(mod_dat$final.est.rich,center=TRUE,scale=TRUE)
mod_dat$fin.est.div<-scale(mod_dat$final.est.div,center=TRUE,scale=TRUE)
mod_dat$turnover<-scale(mod_dat$turnover,center=TRUE,scale=TRUE)

# Richness

rich<-glm(prop.final.est.rich ~ turnover * fin.est.rich * years,
          family='quasibinomial',
          data=mod_dat)
summary(rich)
stepCriterion(rich)

rich.final<-glm(prop.final.est.rich~ fin.est.rich + years + turnover,
                family='quasibinomial',
                data=mod_dat)
rich.sum<-summary(rich.final)
rich.sum
rich.sum$null.deviance
rich.sum$deviance
rich.sum$df.residual
pseudo_r2 <- 1 - (rich.sum$deviance / rich.sum$null.deviance)
pseudo_r2


# Diversity

div<-lm(prop.final.est.div ~ turnover * fin.est.div * years,
        data=mod_dat)
summary(div)
stepCriterion(div)

div.final<-lm(prop.final.est.div~ turnover,
              data=mod_dat)
summary(div.final)

div.null <- lm(prop.final.est.div ~ 1, data = mod_dat)
div.null.dev <- deviance(div.null)
print(div.null.dev)
deviance(div.final)



### Turnover by overlap

signif_model_rich<-glm(signif.rich~ turnover,
                       data=mod_dat,
                       family="binomial")
summary(signif_model_rich)

# no relationship with turnover

signif_model_div<-glm(signif.div~ turnover,
                      data=mod_dat,
                      family="binomial")
summary(signif_model_div)

# prob of overlap decreases with turnover marginally significant
