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
library('lme4')
library("glmtoolbox")
library('lmerTest')

# Load in the formatted clean data
#Make sure the results are correctly configured

full.com<-read.csv(paste0(datapath,'communityResults.csv'))
rich.thresh<-read.csv(paste0(datapath,'richnessThresh.csv'))
div.thresh<-read.csv(paste0(datapath,'diversityThresh.csv'))

# format data for full models

mod_dat<-full.com[which(full.com$year!="full"),] # remove duplication of full data
mod_dat<-mod_dat[which(is.na(mod_dat$turnover)==FALSE),]

mod_dat$years<-scale(mod_dat$years,center=TRUE,scale=TRUE)
mod_dat$fin.est.rich<-scale(mod_dat$final.est.rich,center=TRUE,scale=TRUE)
mod_dat$fin.est.div<-scale(mod_dat$final.est.div,center=TRUE,scale=TRUE)
mod_dat$turnover<-scale(mod_dat$turnover,center=TRUE,scale=TRUE)

mod_dat$obs<-1:nrow(mod_dat)

#### Richness

rich<-glmer(prop.final.est.rich ~ turnover * fin.est.rich * years
            + (1|site) + (1|obs), # random effect of obs included to handle overdispersion, since quasibinomial is not possible in glmer
            family='binomial',
            data=mod_dat,
            control = glmerControl(optimizer = "bobyqa"))
summary(rich)

# random effects explain no variance so they are removed
rich<-glm(prop.final.est.rich ~ turnover * fin.est.rich * years,
            family='quasibinomial',
            data=mod_dat)
summary(rich)
stepCriterion(rich)

rich.final<-glm(prop.final.est.rich ~ years + turnover + fin.est.rich + years:turnover,
          family='quasibinomial',
          data=mod_dat)
summary(rich.final)

# positive effect of the number of years, negative effect of turnover and estimated richness, positive effect of turnoverxyear


#### Diversity

div <- lmer(prop.final.est.div ~turnover * fin.est.div * years
            + (1|site),
            data=mod_dat,
            control = lmerControl(optimizer = "bobyqa"))

summary(div)
step(div)

div.final <- lmer(prop.final.est.div ~turnover + fin.est.div + years + (1 | site) + turnover:years + fin.est.div:years,
            data=mod_dat,
            control = lmerControl(optimizer = "bobyqa"))

summary(div.final)

# negative effect of turnover, positive effect of years, marginally negative effect of div, positive effect of years * turnover and  estimated diversity x years

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


#### Full data only analsis

mod_dat<-full.com[which(full.com$year=="full"),] # fulld data only

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
summary(rich.final)
# Significant decrease with estimated richness, nonsignificant pos trend with years and neg with turnover


# Diversity

div<-lm(prop.final.est.div ~ turnover * fin.est.div * years,
          data=mod_dat)
summary(div)
stepCriterion(div)

div.final<-glm(prop.final.est.div~ turnover,
                data=mod_dat)
summary(div.final)
# negative relationship with turnover

#### reaching thresholds

rich<-lm(y.thresh ~ turnover,
          data=rich.thresh)
summary(rich)
# no relationship between turnover and number of years

div<-lm(y.thresh ~ turnover,
         data=div.thresh)
summary(div)
# number of years increases with turnover

