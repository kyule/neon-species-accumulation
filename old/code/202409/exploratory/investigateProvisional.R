load("/Users/kelsey/Github/neon-species-accumulation/data/NeonData.Robj")

fullData<-NeonData$bet_expertTaxonomistIDProcessed
fullData$year<-year(fullData$collectDate)

sites<-unique(fullData$siteID)
results<-setNames(vector(mode="list",length=length(sites)),sites)



for (i in 1:length(sites)){
  print(sites[i])
  dat<-fullData[which(fullData$siteID==sites[i]),]
  spp<-unique(dat$taxonID)
  years<-unique(dat$year)

  # initiate an incidence data frame
  inc<-data.frame(matrix(ncol=length(years),nrow=length(spp)))
  colnames(inc)<-paste0("year_",years)
  rownames(inc)<-spp
  inc[is.na(inc)]<-0
  
  # Input individuals into incidence data frame
  for (k in 1:nrow(dat)){
    row<-which(rownames(inc)==dat$taxonID[k])
    if (length(row)>0){
      col<-which(colnames(inc)==paste0("year_",dat$year[k]))
      if (length(col)>0){
        inc[row,col]<-1
      }
    }
  }
  
  input<-c(ncol(inc),as.vector(rowSums(inc)))
  richness.out<-iNEXT(input,q=0,datatype='incidence_freq',knot=30,endpoint=ncol(inc)*10)
  results[[i]]<-richness.out
  
  ### Calculate turnover measures
  
  ########### FIX THIS THOUGH, needs to be real individual counts!!
  
  print(paste0("turnover: ",sites[i]))
  turnover.df<-data.frame(dat %>% group_by(year,taxonID) %>% count())#summarise(count=sum(finalIndivCount)))
  #turnover.df<-turnover.df[turnover.df$count>=1,]
  total<-turnover(turnover.df,time.var="year",species.var="taxonID",abundance.var="n",metric="total")
  #appear<-turnover(turnover.df,time.var="year",species.var="taxonID",abundance.var="count",metric="appearance")
  #disappear<-turnover(turnover.df,time.var="year",species.var="taxonID",abundance.var="count",metric="disappearance")
  #turnover<-left_join(total,appear,join_by("year"=="year"))
  #turnover<-left_join(turnover,disappear,join_by("year"=="year"))

  results[[i]]$turnover<-total
  
  # dissimilarity measures
  
  ########### FIX THIS THOUGH, needs to be real individual counts!!
  
  print(paste0("dissimilarity: ",sites[i]))
  dat$siteyear<-paste(dat$siteID,dat$year,sep="_")
  dis.dat<-data.frame(dat %>% group_by(siteyear,taxonID) %>% count())#summarise(n=sum(finalIndivCount)))
  dis.dat<-dis.dat[which(dis.dat$n>0),]
  dis.dat<-acast(dis.dat, siteyear~taxonID, value.var="n")
  dissim<-ess(dis.dat)
  
  results[[i]]$dissim<-mean(dissim)

}

# Pull Estimated Asymptotic and Observed Richness Values out of the results list -- full data only
results<- results[!names(results)=="GUAN"]


asyest.list<-lapply(results,function(item) data.frame(item$AsyEst))
asyest <- bind_rows(asyest.list,.id="site")
richness <- asyest[grep("Richness",rownames(asyest)),]

turnover.list<-lapply(results,function(item) data.frame(item$turnover))
turnover.list <- bind_rows(turnover.list,.id="site")

turnover<-data.frame(turnover.list %>% 
                       group_by(site) %>% 
                       summarise(turnover=mean(total,na.rm=TRUE),years=length(year)))
full.com <- left_join(richness,turnover,join_by("site"=="site"))

full.com$propObs<-full.com$Observed/full.com$Estimator



dissim.list<-lapply(results,function(item) item$dissim)
dissim<-melt(dissim.list)
names(dissim)<-c("dissim","site")
full.com<-left_join(full.com,dissim,join_by("site"=="site"))

summary(lm(Observed~turnover,full.com))
plot(Observed~turnover,full.com)

summary(lm(Estimator~turnover,full.com))
plot(Estimator~turnover,full.com)

summary(glm(propObs~turnover,full.com,family=binomial(link="logit")))
ggplot(full.com, aes(x = turnover, y = propObs, color=as.numeric(Estimator))) +
  geom_point(aes(size=years)) + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), color = "black") +  # Add linear regression line
  labs(x = "Mean Species Turnover", y = "Prop. Estimated Species Observed") +
  theme_minimal() +
  scale_color_viridis_c(option = "D",name="Est. Richness")
# negative relationship between turnover and proportion of estimated species richness we have observed

summary(glm(propObs~dissim,full.com,family=binomial(link="logit")))
ggplot(full.com, aes(x = dissim, y = propObs, color=as.numeric(Estimator))) +
  geom_point(aes(size=years)) + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), color = "black") +  # Add linear regression line
  labs(x = "Mean Dissimilarity", y = "Prop. Estimated Species Observed") +
  theme_minimal() +
  scale_color_viridis_c(option = "D",name="Est. Richness")

# Take into account the number of traps that have been sampled

trap.sum <- data.frame(completeness %>% group_by(siteID) %>% summarise(traps=sum(traps)))
full.com<-left_join(full.com,trap.sum,join_by("site"=="siteID"))

# Start some basic analyses
summary(lm(Observed~turnover,full.com))
plot(Observed~turnover,full.com)
# No relationship between observed and turnover

summary(lm(Estimator~turnover,full.com))
plot(Estimator~turnover,full.com)

# Estimated increases with turnover

summary(glm(propObs~turnover,full.com,family=binomial(link="logit")))

ggplot(full.com, aes(x = turnover, y = propObs, color=as.numeric(Estimator))) +
  geom_point(aes(size=years)) + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), color = "black") +  # Add linear regression line
  labs(x = "Mean Species Turnover", y = "Prop. Estimated Species Observed") +
  theme_minimal() +
  scale_color_viridis_c(option = "D",name="Est. Richness")
# negative relationship between turnover and proportion of estimated species richness we have observed


summary(glm(propObs~Estimator,full.com,family=binomial(link="logit")))
# negative relationship between estimated spp richness and proportion of estimated species richness we have observed
ggplot(full.com, aes(x = Estimator, y = propObs, color=as.numeric(turnover))) +
  geom_point(size=3) + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), color = "black") +  # Add linear regression line
  labs(y = "Prop. Estimated Species Observed", x = "Estimated Richness") +
  theme_minimal() +
  scale_color_viridis_c(option = "D",name="turnover")
# negative relationship between the proportion observed and the estimated richness overall

summary(glm(propObs~years+turnover,full.com,family=binomial(link="logit")))
#There is no relationship between the proportion observed and the number of years sampled overall

summary(lm(Estimator~years,full.com))
#estimated richness increases with years of sampling

fullmodel<-glm(propObs~years*Estimator*turnover,full.com,family=binomial(link="logit"))
stepaic<-stepAIC(fullmodel,direction="both")
summary(stepaic)


### Find sampling required to reach 90% diversity

# Pull Estimated and Observed Richness by sampling Values out of the results list -- full community only

inext.list<-lapply(results,function(item) data.frame(item$iNextEst$size_based))
inext <- bind_rows(inext.list,.id="site")
inext <- inext %>% filter(Order.q==0)

thresh90<-data.frame(site=full.com$site,thresh=full.com$Estimator*0.90)


# Find intersections
find_intersections <- function(df,thresh,metric) {
  site_data <- df[df$site == thresh$site, ]
  
  # Loop through data and check for crossings
  for (i in 1:(nrow(site_data) - 1)) {
    if ((site_data$qD[i] < thresh$thresh && site_data$qD[i+1] > thresh$thresh) ||
        (site_data$qD[i] > thresh$thresh && site_data$qD[i+1] < thresh$thresh)) {
      
      # Linear interpolation to find the exact crossing point
      if (metric=="t"){
        t1 <- site_data$t[i]
        t2 <- site_data$t[i+1]
        qD1 <- site_data$qD[i]
        qD2 <- site_data$qD[i+1]
        crossing_t <- t1 + (thresh$thresh - qD1) * (t2 - t1) / (qD2 - qD1)
        
        return(data.frame(site = thresh$site, t = crossing_t, qD = thresh$thresh))
      }
      
      if (metric=="byYears")
        y1 <- site_data$byYears[i]
      y2 <- site_data$byYears[i+1]
      qD1 <- site_data$qD[i]
      qD2 <- site_data$qD[i+1]
      crossing_t <- y1 + (thresh$thresh - qD1) * (y2 - y1) / (qD2 - qD1)
      
      return(data.frame(site = thresh$site, t = crossing_t, qD = thresh$thresh))
    }
  }
  return(NULL)
}

# Get all intersections 
intersections <- do.call(rbind, lapply(1:nrow(thresh90), function(i) {
  find_intersections(inext, thresh90[i, ],"t")
}))

full.com<-left_join(full.com,intersections,join_by("site"=="site"))

ggplot(inext, aes(x = t, y = qD, color=site , group=site)) +
  geom_line() +
  labs(x = "Number of Years", y = "Estimated Richness") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1") +
  geom_hline(data = thresh90, aes(yintercept = thresh, color = site), linetype = "dashed") +
  geom_point(data = intersections, aes(x = t, y = qD), color = "darkgrey", size = 2) 

fullmodel<-lm(t~Estimator*turnover,full.com)
stepaic<-stepAIC(fullmodel,direction="both")
summary(stepaic)

# It takes longer the more species are estimated to be there, #
# but that fact is  dissimilar communities with more 
