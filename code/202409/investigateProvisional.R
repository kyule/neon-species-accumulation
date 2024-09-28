load("/Users/kelsey/Github/neon-species-accumulation/data/testIncludingProvisional/NeonData.Robj")

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



