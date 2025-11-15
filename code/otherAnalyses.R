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
library("patchwork")
library('viridis')
library('ggplot2')

## Plot proportion observed as a function of turnover: full data only


rich.full <- ggplot(full.com[which(full.com$year == "full"),], aes(x = turnover, y = prop.final.est.rich, color = as.numeric(Estimator.rich))) +
  geom_point(aes(size = years)) +
  geom_smooth(method = "glm", method.args = list(family = "quasibinomial"), color = "black") +
  labs(x = "", y = "Observed/Estimated Richness") +
  theme_minimal() +
  annotate("text", x = -Inf, y = Inf, label = "a", fontface = "bold", hjust = -0.2, vjust = 1.3, size = 6) +
  scale_color_viridis_c(option = "D", name = "Est. richness",limits = rich_limits)

div.full <- ggplot(full.com[which(full.com$year == "full"),], aes(x = turnover, y = prop.final.est.div, color = as.numeric(Estimator.div))) +
  geom_point(aes(size = years)) +
  geom_smooth(method = "lm", color = "black") +
  labs(x = "Mean Species Turnover", y = "Observed/Estimated Diversity") +
  theme_minimal() +
  annotate("text", x = -Inf, y = Inf, label = "b", fontface = "bold", hjust = -0.2, vjust = 1, size = 6) +
  scale_color_viridis_c(option = "D", name = "Est. diversity", limits = div_limits) +
  guides(size = "none") 

combined.full <- rich.full / div.full + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "right")

print(combined.full)


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

# format data for full models

mod_dat<-full.com[which(full.com$year!="full"),] # remove duplication of full data
mod_dat<-mod_dat[which(is.na(mod_dat$turnover)==FALSE),]

mod_dat$years<-scale(mod_dat$years,center=TRUE,scale=TRUE)
mod_dat$fin.est.rich<-scale(mod_dat$final.est.rich,center=TRUE,scale=TRUE)
mod_dat$fin.est.div<-scale(mod_dat$final.est.div,center=TRUE,scale=TRUE)
mod_dat$turnover<-scale(mod_dat$turnover,center=TRUE,scale=TRUE)

mod_dat$obs<-1:nrow(mod_dat)

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



### Violin plots

# Violin plots of turnover vs. overlap with estimator +/- se: all data only


richplot<-ggplot(full.com[which((full.com$year != "full") & !is.na(full.com$turnover)),],aes(x=turnover, y=as.factor(signifText.rich)))+
  labs(y = "Richness", x = "Mean Species Turnover") +
  geom_violin() + stat_summary(
    fun = "mean", geom = "point", shape = 20, size = 3, color = "black")+
  theme_minimal()   + theme(
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(),    
    axis.title.y = element_text(size = 18),  
    axis.text.y = element_text(size = 16)     
  )

divplot<-ggplot(full.com[which((full.com$year != "full") & !is.na(full.com$turnover)),],aes(x=turnover, y=as.factor(signifText.div)))+
  labs(y = "Diversity", x = "Mean Species Turnover") +
  geom_violin() + stat_summary(
    fun = "mean", geom = "point", shape = 20, size = 3, color = "black")+
  theme_minimal()   + theme(
    axis.title = element_text(size = 18),  
    axis.text = element_text(size = 16)  
  )

combined_plot <- richplot/divplot
combined_plot

# Violin plots of turnover vs. overlap with estimator +/- se: full data only


richplot<-ggplot(full.com[which(full.com$year == "full"),],aes(x=turnover, y=as.factor(signifText.rich)))+
  labs(y = "Richness", x = "Mean Species Turnover") +
  geom_violin() + stat_summary(
    fun = "mean", geom = "point", shape = 20, size = 3, color = "black")+
  theme_minimal()   + theme(
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(),    
    axis.title.y = element_text(size = 18),  
    axis.text.y = element_text(size = 16)     
  )

divplot<-ggplot(full.com[which(full.com$year == "full"),],aes(x=turnover, y=as.factor(signifText.div)))+
  labs(y = "Diversity", x = "Mean Species Turnover") +
  geom_violin() + stat_summary(
    fun = "mean", geom = "point", shape = 20, size = 3, color = "black")+
  theme_minimal()   + theme(
    axis.title = element_text(size = 18),  
    axis.text = element_text(size = 16)  
  )

combined_plot <- richplot/divplot
combined_plot

# Violin plots of turnover vs. overlap with estimator +/- se: all data 


richplot<-ggplot(full.com[which((full.com$year != "full") & !is.na(full.com$turnover)),],aes(x=turnover, y=as.factor(signifText.rich)))+
  labs(y = "Richness", x = "Mean Species Turnover") +
  geom_violin() + stat_summary(
    fun = "mean", geom = "point", shape = 20, size = 3, color = "black")+
  theme_minimal()   + theme(
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(),    
    axis.title.y = element_text(size = 18),  
    axis.text.y = element_text(size = 16)     
  )

divplot<-ggplot(full.com[which((full.com$year != "full") & !is.na(full.com$turnover)),],aes(x=turnover, y=as.factor(signifText.div)))+
  labs(y = "Diversity", x = "Mean Species Turnover") +
  geom_violin() + stat_summary(
    fun = "mean", geom = "point", shape = 20, size = 3, color = "black")+
  theme_minimal()   + theme(
    axis.title = element_text(size = 18),  
    axis.text = element_text(size = 16)  
  )

combined_plot <- richplot/divplot
combined_plot

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


####### Turnover figures
# Load in the formatted clean data
#Make sure the results are correctly configured

load(file=paste0(datapath,"iNEXTandTurnoverResults.Robj"))
full.com<-read.csv(paste0(datapath,'communityResults.csv'))
full.com<-full.com[which(full.com$year=="full"),]

# Remove GUAN due to low data availability
results<- results[!names(results)=="GUAN"]

# Grab turnover

turnover.list<-lapply(results,function(item) data.frame(item$turnover))
turnover.list <- bind_rows(turnover.list,.id="site")
turnover<-data.frame(turnover.list %>% 
                       group_by(site) %>% 
                       summarise(turnover=mean(total,na.rm=TRUE),years=length(year)))


# Pull Estimated and Observed Richness and diversity by sampling Values out of the results list -- full community only

inext.list<-lapply(results,function(item) data.frame(item$full$out$iNextEst$size_based))
inext <- bind_rows(inext.list,.id="site")
inext.rich <- inext %>% filter(Order.q==0)
inext.div <- inext %>% filter(Order.q==1)

thresh90.rich<-data.frame(site=full.com$site,thresh=full.com$Estimator.rich*0.90)
thresh90.div<-data.frame(site=full.com$site,thresh=full.com$Estimator.div*0.90)
full.com$trapAvg<-full.com$traps/full.com$years
thresh90.rich<-left_join(thresh90.rich,full.com[,c("site","trapAvg")],join_by("site"=="site"))
thresh90.div<-left_join(thresh90.div,full.com[,c("site","trapAvg")],join_by("site"=="site"))


## Find where the thresholds is crossed and scale years in inext data frame
thresh90.rich$t.thresh<-NA
inext.rich$y<-NA

for (i in 1:nrow(thresh90.rich)){
  site_data<-inext.rich[which(inext.rich$site==thresh90.rich$site[i]),]
  inext.rich$y[which(inext.rich$site==thresh90.rich$site[i])]<-inext.rich$t[which(inext.rich$site==thresh90.rich$site[i])]/thresh90.rich$trapAvg[i]
  diff<-site_data$qD-thresh90.rich$thresh[i]
  
  for (j in 1:length(diff)){
    if (diff[j]>0){break}
  }
  
  v<-site_data[c(j-1,j),c("t","qD")]
  thresh90.rich$t.thresh[i]<-v$t[1] + (thresh90.rich$thresh[i] - v$qD[1]) * (v$t[2] - v$t[1]) / (v$qD[2] - v$qD[1])
  
}

thresh90.div$t.thresh<-NA
inext.div$y<-NA

for (i in 1:nrow(thresh90.div)){
  site_data<-inext.div[which(inext.div$site==thresh90.div$site[i]),]
  inext.div$y[which(inext.div$site==thresh90.div$site[i])]<-inext.div$t[which(inext.div$site==thresh90.div$site[i])]/thresh90.div$trapAvg[i]
  diff<-site_data$qD-thresh90.div$thresh[i]
  
  for (j in 1:length(diff)){
    if (diff[j]>0){break}
  }
  
  v<-site_data[c(j-1,j),c("t","qD")]
  thresh90.div$t.thresh[i]<-v$t[1] + (thresh90.div$thresh[i] - v$qD[1]) * (v$t[2] - v$t[1]) / (v$qD[2] - v$qD[1])
  
}


# Calculate the thresholds based on year instead 
thresh90.rich$y.thresh<-thresh90.rich$t.thresh/thresh90.rich$trapAvg
thresh90.div$y.thresh<-thresh90.div$t.thresh/thresh90.div$trapAvg

full.com.rich<-left_join(full.com,thresh90.rich,join_by("site"=="site"))
full.com.div<-left_join(full.com,thresh90.div,join_by("site"=="site"))

# Bring in field site data, join other data for easy plotting

sites<-read.csv(paste0(datapath,"NEON_Field_Site_Metadata_20240926.csv"))
domains<-data.frame(ID=c('D01',	'D02',	'D03',	'D04',	'D05',	'D06',	'D07',	'D08',	'D09',	'D10',	'D11',	'D12',	'D13',	'D14',	'D15',	'D16',	'D17',	'D18',	'D19',	'D20'),
                    domainName=c('Northeast',	'Mid-Atlantic',	'Southeast',	'Atlantic Neotropical',	'Great Lakes',	'Prairie Peninsula',	'Central Plains',	'Ozarks Complex',	'Northern Plains',	'Southern Plains',	'Southern Rockies/Colorado Plateau',	'Desert Southwest',	'Pacific Southwest',	'Northern Rockies',	'Great Basin',	'Pacific Northwest',	'California',	'Tundra',	'Taiga',	'Pacific Tropical'))
sites<-left_join(sites,domains,join_by("field_domain_id"=="ID"))
inext<-left_join(inext,sites,join_by("site"=="field_site_id"))
thresh90.rich<-left_join(thresh90.rich,sites,join_by("site"=="field_site_id"))
thresh90.div<-left_join(thresh90.div,sites,join_by("site"=="field_site_id"))
thresh90.rich<-left_join(thresh90.rich,turnover,join_by("site"=="site"))
thresh90.div<-left_join(thresh90.div,turnover,join_by("site"=="site"))

inext.rich<-left_join(inext.rich,thresh90.rich,join_by("site"=="site"))
inext.div<-left_join(inext.div,thresh90.div,join_by("site"=="site"))

# Find x and ylims
xlim.rich<-ceiling(max(thresh90.rich$y.thresh))
ylim.rich<-ceiling(max(full.com$Estimator.rich))
xlim.div<-ceiling(max(thresh90.div$y.thresh))
ylim.div<-ceiling(max(full.com$Estimator.div))

#Plot the accumulation curves and estimated number of years
ggplot(inext.rich, aes(x = y, y = qD,group=site,color=as.numeric(turnover))) +
  geom_line(size=1) +
  facet_wrap(~ domainName) +
  labs(x = "Years", y = "Estimated Richness") +
  theme_minimal() +
  geom_hline(data = thresh90.rich, aes(yintercept = thresh), linetype = "dashed") +
  geom_point(data = thresh90.rich, aes(x = y.thresh, y = thresh), color = "darkgrey", size = 2) +
  geom_point(data = inext.rich[which(inext.rich$Method=="Observed"),],aes(x= y, y=qD),color="black", size = 2) +
  ylim(0,ylim.rich) +
  xlim(0,xlim.rich) +
  scale_color_viridis_c(option = "D",name="Avg. turnover")

ggplot(inext.div, aes(x = y, y = qD,group=site,color=as.numeric(turnover))) +
  geom_line(size=1) +
  facet_wrap(~ domainName) +
  labs(x = "Years", y = "Estimated Diversity") +
  theme_minimal() +
  geom_hline(data = thresh90.div, aes(yintercept = thresh), linetype = "dashed") +
  geom_point(data = thresh90.div, aes(x = y.thresh, y = thresh), color = "darkgrey", size = 2) +
  geom_point(data = inext.div[which(inext.div$Method=="Observed"),],aes(x= y, y=qD),color="black", size = 2) +
  ylim(0,ylim.div) +
  xlim(0,xlim.div) +
  scale_color_viridis_c(option = "D",name="Avg. turnover")
