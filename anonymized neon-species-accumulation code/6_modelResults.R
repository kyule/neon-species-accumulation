### This is the primary file for analyzing NEON species accumulation data

#### BEFORE RUNNING
# Users must define their own paths

# datapath<-"user defined path"
# codepath<-"user defined path"



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

rich.norand<-glmer(prop.final.est.rich ~ turnover * fin.est.rich * years
            + (1|obs), # random effect of obs included to handle overdispersion, since quasibinomial is not possible in glmer
            family='binomial',
            data=mod_dat,
            control = glmerControl(optimizer = "bobyqa"))
summary(rich.norand)


rich<-glm(prop.final.est.rich ~ turnover * fin.est.rich * years,
            family='quasibinomial',
            data=mod_dat)
summary(rich)
stepCriterion(rich)

rich.final<-glm(prop.final.est.rich ~ years + turnover + fin.est.rich + years:turnover,
          family='quasibinomial',
          data=mod_dat)
summary(rich.final)


#### Diversity

div <- lmer(prop.final.est.div ~turnover * fin.est.div * years
            + (1|site),
            data=mod_dat,
            control = lmerControl(optimizer = "bobyqa"))

div.sum<-summary(div)

div.norand <- lm(prop.final.est.div ~turnover * fin.est.div * years,
            data=mod_dat)

div.norand.sum<-summary(div.norand)
logLik(div)
logLik(div.norand)

step(div.norand)

div.final <- lm(prop.final.est.div ~ turnover + fin.est.div + years + 
                  turnover:fin.est.div + turnover:years + fin.est.div:years,
            data=mod_dat)
summary(div.final)

div.null <- lm(prop.final.est.div ~ 1, data = mod_dat)
div.null.dev <- deviance(div.null)
print(div.null.dev)
deviance(div.final)




