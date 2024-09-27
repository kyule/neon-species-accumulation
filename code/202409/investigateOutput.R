which(lapply(results, function(item) item$out)=="Error")

results$DSNY$inc_freq
results$DSNY$inc_freq

turnover.results<-lapply(results,function(item) item$turnover$total)

turnover.df <- data.frame(
  name = rep(names(turnover.results), times = sapply(turnover.results, length)),
  value = unlist(turnover.results)
)
turn<-turnover.df %>% group_by(name) %>% summarise(turn=mean(value))

inext.results<-lapply(results,function(item) item$out)
inext.results<-lapply(inext.results[1:5],function(item) item$AsyEst)
inext.results<-lapply(inext.results,function(item) item[which(item$Assemblage=="full" & item$Diversity=="Species richness"),])
inext.results<-lapply(inext.results,function(item) item[,3]/item[,4])

inext.df <- data.frame(
  name = rep(names(inext.results), times = sapply(inext.results, length)),
  value = unlist(inext.results)
)
full_join(inext.df,turn,join_by("name"=="name"))[1:5,]->examp

plot(turn~value,examp)

inext.results<-lapply(results,function(item) item$out)
inext.results<-lapply(inext.results[1:5],function(item) item$AsyEst)
inext.results<-lapply(inext.results,function(item) item[which(item$Assemblage=="full" & item$Diversity=="Species richness"),])
inext.results<-lapply(inext.results,function(item) item[,4])

inext.df <- data.frame(
  name = rep(names(inext.results), times = sapply(inext.results, length)),
  value = unlist(inext.results)
)
full_join(inext.df,turn,join_by("name"=="name"))[1:5,]->examp

plot(turn~value,examp)



ggiNEXT(iNEXT(results$BART[[5]]$inc_freq,q=c(0),datatype='incidence_freq',knot=20,endpoint=ncol(results$BART[[5]]$inc)*3))
