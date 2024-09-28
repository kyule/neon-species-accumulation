### Find sampling required to reach 90% diversity

# Pull Estimated and Observed Richness by sampling Values out of the results list -- full community only

inext.list<-lapply(results,function(item) data.frame(item$full$out$iNextEst$size_based))
inext <- bind_rows(inext.list,.id="site")
inext <- inext %>% filter(Order.q==0)

thresh90<-data.frame(site=full.com$site,thresh=full.com$Estimator*0.90)
full.com$trapAvg<-full.com$traps/full.com$years
thresh90<-left_join(thresh90,full.com[,c("site","trapAvg")],join_by("site"=="site"))

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




ggplot(inext, aes(x = t, y = qD, color=site , group=site)) +
  geom_line() +
  labs(x = "Number of Traps", y = "Estimated Richness") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1") +
  geom_hline(data = thresh90, aes(yintercept = thresh, color = site), linetype = "dashed") +
  geom_point(data = intersections, aes(x = t, y = qD), color = "darkgrey", size = 2) 

#Calculate the number of years based on recent per year sampling effort

thresh90<-left_join(thresh90,intersections[c("site","t")],join_by("site"=="site"))
thresh90$yearsReq<-thresh90$t/thresh90$trapno
thresh90

# Scale the X axis by years of sampling

inext$byYears<-NA
for (i in 1:nrow(thresh90)){
  inext$byYears[which(inext$site==thresh90$site[i])]<-inext$t[which(inext$site==thresh90$site[i])]/thresh90$trapno[i]
}

intersections <- do.call(rbind, lapply(1:nrow(thresh90), function(i) {
  find_intersections(inext, thresh90[i, ],"byYears")
}))


ggplot(inext, aes(x = byYears, y = qD, color=site , group=site)) +
  geom_line() +
  labs(x = "Number of Years", y = "Estimated Richness") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1") +
  geom_hline(data = thresh90, aes(yintercept = thresh, color = site), linetype = "dashed") +
  geom_point(data = intersections, aes(x = t, y = qD), color = "black", size = 2) 



# Turnover relates to number of years?

thresh90<-left_join(thresh90,turnover,join_by("site"=="site"))

summary(lm(yearsReq~turnover,thresh90))






