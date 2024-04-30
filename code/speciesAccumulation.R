### This is the primary file for analyzing NEON species accumulation data

#### BEFORE RUNNING
# Users must define their own paths

# datapath<-"user defined path"
# codepath<-"user defined path"

#And users must have configured the DownloadNEONData.R file as desired for their work

# And users should remove below line that loads in personal paths
source("/Users/kelsey/Github/neon-species-accumulation/configini.R")

# Load some necessary packages
library("ggplot2")
library("dplyr")

# Load in the formatted data tables and sampling effort

source(paste0(codepath,"DataCleaning.R"))
source(paste0(codepath,"SampleEffort.R"))


# Join expert and para tables

betIDs<-left_join(para,expert,join_by("individualID"=="individualID"),suffix=c(".para",".expert"))
betIDs[which(betIDs$uid.para %in% betIDs$uid.para[which(duplicated(betIDs$uid.para))]),]
# Duplicates are potentially re-identifications by different experts of the same beetles, leaving all data to be conservative, does not affect overall presence of species at any site

#subset to when expert exists
betIDsExp<-betIDs[which(is.na(betIDs$scientificName.expert)==FALSE),]
mismatch<-betIDsExp[which(betIDsExp$scientificName.para!=betIDsExp$scientificName.expert),]
#see that many mismatches exist, so need to propagate to rest of IDs


# Put expert ID where exist into sorting table

for (i in 1:nrow(expert)){
  #what sorting table beetles are associated with this pinned beetle identification
  sub<-betIDs$subsampleID[which(betIDs$individualID==expert$individualID[i])]
  sp<-betIDs$taxonID.para[which(betIDs$individualID==expert$individualID[i])]
  sortRecs<-which(sort$subsampleID==sub & sort$taxonID==sp)
  sort$taxonID[i]<-expert$taxonID[i]
}

# subset sort table for ease of use
sort<-sort[,which(names(sort) %in% c("uid.x","subsampleID","sampleID","taxonID","individualCount"))]

# join sort table with taxonomy table for further analysis
taxa<-taxa[which(names(taxa) %in% c("taxonID","scientificName","taxonRank","family","subfamily","tribe","genus","subgenus","specificEpithet","infraspecificEpithet"))]
sort<-left_join(sort,taxa,join_by("taxonID"=="taxonID"),suffix=c("sort","taxa"))

#remove non-carabids, replace OTHE with CARSP14
sort<-sort[which(sort$family=="Carabidae"),]
sort$taxonID[which(sort$taxonID=="OTHE")]<-"CARSP14"
sort$scientificName[which(sort$taxonID=="CARSP14")]<-"Carabidae sp."

# investigate taxonomic ranks in data
sort %>% group_by(taxonRank) %>% count()
# about 86.1k out of 89.9k identified to species or subspecies (>95%), 88.6k out of 89.9k identified to at least genus (98.5%)

# Want to figure out what to determine a unique scientific name. 
####Let's ignore subspecies and give those with a subspecies ID the taxon at the species level

sort$analysisSciName<-sort$scientificName 
subsp<-sort[which(sort$taxonRank=="subspecies"),]
subsp$analysisSciName<-paste(sapply(strsplit(subsp$scientificName," "),"[",1),sapply(strsplit(subsp$scientificName," "),"[",2),sep=" ")

for (i in 1:nrow(subsp)){
  a<-which(sort$uid.x==subsp$uid.x[i])
  sort$analysisSciName[a]<-subsp$analysisSciName[i]
}

# decided to leave other taxonomic rank issues for now
#subset field table to values of interest and combine with sort table
field<-field[,which(names(field) %in% c("siteID","year","nlcdClass","sampleID"))]
data<-left_join(sort,field,join_by("sampleID"=="sampleID"))
# all data have matches to field data :)

# investigate taxonRank vs year
rankyears<-data %>% group_by(year,taxonRank) %>% count()
years<-data %>% group_by(year) %>% count()
rankyearprop<-left_join(rankyears,years,join_by("year"=="year"))
rankyearprop$prop<-rankyearprop$n.x/rankyearprop$n.y

ggplot(data=rankyearprop[-which(rankyearprop$taxonRank %in% c("species","subspecies")),], aes(x=as.numeric(year), y=prop, group=taxonRank)) +
  geom_line(aes(color=taxonRank))

# individuals being identified only to family drops rapidly over time. Subgenus has a strong peak in 2017
# Decide to throw out individuals identified only to family
data<-data[-which(data$taxonRank=="family"),]

# plot accumulation curve

data$year<-as.numeric(data$year)

cumulativeSpp.Year<- function(site){
  dat.site<-data[which(data$siteID==site),]
  dat.site$year<-dat.site$year-min(dat.site$year)+1
  years<-min(dat.site$year):max(dat.site$year)
  nlcds<-unique(dat.site$nlcdClass)
  cumul<-data.frame(year=c(),nlcd=c(),spp=c())
  for (i in 1:length(nlcds)){
    dat.nlcd<-dat.site[which(dat.site$nlcdClass==nlcds[i]),]
    for (j in 1:length(years)){
      acum<-length(unique(dat.nlcd$analysisSciName[which(dat.nlcd$year<=years[j])]))
      cumul<-rbind(cumul,data.frame(year=years[j],nlcd=nlcds[i],spp=acum))
    }
  }
  return(cbind(cumul,data.frame(siteID=site)))

  return(print(plotSite))
}

HARV<-cumulativeSpp.Year("HARV")
plotHARV<-ggplot(HARV, aes(x=year, y=spp, group=nlcd)) +
  geom_line(aes(color=nlcd))+
  geom_point(aes(color=nlcd)) + theme_classic()
plotHARV

SRER<-cumulativeSpp.Year("SRER")
plotSRER<-ggplot(SRER, aes(x=year, y=spp, group=nlcd)) +
  geom_line(aes(color=nlcd))+
  geom_point(aes(color=nlcd)) + theme_classic()
plotSRER

GUAN<-cumulativeSpp.Year("GUAN")
plotGUAN<-ggplot(GUAN, aes(x=year, y=spp, group=nlcd)) +
  geom_line(aes(color=nlcd))+
  geom_point(aes(color=nlcd)) + theme_classic()
plotGUAN

LAJA<-cumulativeSpp.Year("LAJA")
plotLAJA<-ggplot(LAJA, aes(x=year, y=spp, group=nlcd)) +
  geom_line(aes(color=nlcd))+
  geom_point(aes(color=nlcd)) + theme_classic()
plotLAJA

#model <- drm(SRER$spp ~ SRER$year, type="Poisson",fct = DRC.asymReg())





### need to read this: https://bookdown.org/mike/data_analysis/non-linear-regression.html





#### current working stopping point



# subset sort to important variables and join with betIDs

sort<-sort[,which(names(sort) %in% c("subsampleID","sampleID","taxonID","scientificName","taxonRank","morphospeciesID","individualCount"))]
names(sort)[which(names(sort)=="taxonID")]<-"sort_taxonID"
names(sort)[which(names(sort)=="scientificName")]<-"sort_scientificName"
names(sort)[which(names(sort)=="taxonRank")]<-"sort_taxonRank"
names(sort)[which(names(sort)=="morphospeciesID")]<-"sort_morphospeciesID"

betsortIDs<-full_join(betIDs,sort,join_by("subsampleID"=="subsampleID"))

### subset field to important variables and join with betsortIDs

field<-field[which(field$sampleCollected=="Y"),]
field$siteHabitat<-paste(field$siteID,field$nlcdClass,sep="_")
field<-field[,which(names(field) %in% c("sampleID","siteID","collectDate","siteHabitat","trappingDays"))]

data<-full_join(field,betsortIDs,join_by("sampleID"=="sampleID"))

#Make a table of trapping days x siteHabitat x year
data$year<-format(as.Date(data$collectDate),'%Y')
summarySamplingEffort<-data[which(duplicated(data$sampleID)==FALSE),]

summarySamplingEffort










for (site in allsites){
  print(site)
  sitedata<-NeonData$bet_expertTaxonomistIDProcessed[NeonData$bet_expertTaxonomistIDProcessed$siteID == site,]
  to2014<-length(unique(subset(sitedata, format(as.Date(collectDate), "%Y")<=2014)$taxonID))
  to2015<-length(unique(subset(sitedata, format(as.Date(collectDate), "%Y")<=2015)$taxonID))
  to2016<-length(unique(subset(sitedata, format(as.Date(collectDate), "%Y")<=2016)$taxonID))
  to2017<-length(unique(subset(sitedata, format(as.Date(collectDate), "%Y")<=2017)$taxonID))
  to2018<-length(unique(subset(sitedata, format(as.Date(collectDate), "%Y")<=2018)$taxonID))
  to2019<-length(unique(subset(sitedata, format(as.Date(collectDate), "%Y")<=2019)$taxonID))
  to2020<-length(unique(subset(sitedata, format(as.Date(collectDate), "%Y")<=2020)$taxonID))
  to2021<-length(unique(subset(sitedata, format(as.Date(collectDate), "%Y")<=2021)$taxonID))
  to2022<-length(unique(subset(sitedata, format(as.Date(collectDate), "%Y")<=2022)$taxonID))
  
  newrow<-c(site,to2014,to2015,to2016,to2017,to2018,to2019,to2020,to2021,to2022)
  output[nrow(output) +1,]<-newrow
}

#transpose the data frame to swap rows and columns and rename the columns
plot_data<-data.frame(t(output[-1]))
colnames(plot_data)<-output[, 1]
plot_data$year<-c(2014,2015,2016,2017,2018,2019,2020,2021,2022)
View(plot_data)

plot(plot_data$year,plot_data$SRER)
plot(plot_data$year,plot_data$JORN)


ggplot(data=plot_data, aes(x=as.integer(year), y=as.integer(SRER)) )+
  geom_line()

ggplot(data=plot_data, aes(x=year)) + 
  geom_line(aes(y=as.integer(SRER)), color="red") +
  geom_line(aes(y=as.integer(JORN)), color="blue")


install.packages("reshape2")
library("reshape2")
plotdata2<-melt(plot_data , id.vars="year", variable.name="plot")
plotdata3<-plotdata2[plotdata2$value != 0,]
fullplot<-ggplot(data=plotdata3, aes(x=year, y=as.integer(value))) +
  geom_line(aes(colour = plot)) +
  ggtitle("Cumulative Species Diversity by Year") +
  labs(x = "Year", y = "Cumulative Species") +
  theme_bw() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))
ggsave("C:/Users/avacl/OneDrive/Desktop/fullplot.png",fullplot,width=8.09,height=5)


plotdata3$group<-NA
plotdata3$group[plotdata3$plot %in% c("BART", "HARV", "SCBI", "DSNY", "JERC", "OSBS", "UNDE", "ORNL", "TALL", "WOOD", "CPER", "STER", "ONAQ")]<-1
plotdata3$group[plotdata3$plot %in% c("UNDE", "ORNL", "TALL", "WOOD", "CPER", "STER", "ONAQ")]<-2
View(plotdata3)

plotdata3$group[plotdata3$plot %in% c("BLAN", "SERC", "STEI", "TREE", "KONZ", "UKFS", "GRSM", "DELA", "OAES", "MOAB", "NIWO", "JORN", "HEAL")]<-3
plotdata3$group[plotdata3$plot %in% c("GRSM", "DELA", "OAES", "MOAB", "NIWO", "JORN", "HEAL")]<-4

plotdata3$group[plotdata3$plot %in% c("GUAN", "LAJA", "LENO", "NOGP", "CLBJ", "SRER", "ABBY", "SJER", "DEJU")]<-5
plotdata3$group[plotdata3$plot %in% c("SRER", "ABBY", "SJER", "DEJU")]<-6

plotdata3$group[plotdata3$plot %in% c("KONA", "MLBS", "DCFS", "RMNP", "BARR", "TOOL", "BONA")]<-7

plotdata3$group[plotdata3$plot %in% c("YELL", "WREF", "SOAP", "TEAK", "PUUM")]<-8

g1<-ggplot(data=subset(plotdata3, group==1 & !is.na(plot)), aes(x=year, y=as.integer(value))) +
  geom_line(aes(colour = plot), linewidth = 1.5) +
  scale_color_viridis(option="H", discrete = TRUE) +
  #ggtitle("Cumulative Species Diversity by Year") +
  labs(x = "Year", y = "Cumulative Species") +
  theme_bw() +
  theme(legend.position="right") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("C:/Users/avacl/OneDrive/Desktop/group1.png",g1,width=8.09,height=5)
#repeat for all other graphs

g2<-ggplot(data=subset(plotdata3, group==2 & !is.na(plot)), aes(x=year, y=as.integer(value))) +
  geom_line(aes(colour = plot), linewidth = 1.5) +
  scale_color_viridis(option="H", discrete = TRUE) +
  #ggtitle("Cumulative Species Diversity by Year") +
  labs(x = "Year", y = "Cumulative Species") +
  theme_bw() +
  theme(legend.position="right") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("C:/Users/avacl/OneDrive/Desktop/group2.png",g2,width=8.09,height=5)

#example code from email
install.packages("dplyr")
library("dplyr")
a<-NeonData$bet_sorting[NeonData$bet_sorting$siteID == 'BART',]
#remove invert bycatch, herps, etc
b<-a[a$sampleType == 'carabid',]

c<-b %>%
  group_by(year = format(as.Date(setDate),"%Y")) %>%
  summarise(totalCarabids = sum(individualCount), trapsPresent = length(subsampleCode))
View(c)

d<- b %>%
  group_by(scientificName, year = format(as.Date(setDate),"%Y")) %>%
  summarise(totalCarabids = sum(individualCount), trapsPresent = length(subsampleCode))
View(d)



#testing
a<-NeonData$bet_sorting[NeonData$bet_sorting$siteID == 'SRER',]
#remove invert bycatch, herps, etc
b<-a[a$sampleType == 'carabid',]

c<-b %>%
  group_by(year = format(as.Date(setDate),"%Y")) %>%
  summarise(totalCarabids = sum(individualCount), trapsPresent = length(subsampleCode))
View(c)

d<- b %>%
  group_by(scientificName, year = format(as.Date(setDate),"%Y")) %>%
  summarise(totalCarabids = sum(individualCount), trapsPresent = length(subsampleCode))
View(d)

sp2016<-unique(d$scientificName[d$year == 2016])
sp2017<-unique(d$scientificName[d$year == 2017])
sp2018<-unique(d$scientificName[d$year == 2018])
sp2019<-unique(d$scientificName[d$year == 2019])
sp2021<-unique(d$scientificName[d$year == 2021])

new2017<-subset(d, year == 2017 & !(scientificName %in% sp2016))
new2018<-subset(d, year == 2018 & !(scientificName %in% sp2017))
new2019<-subset(d, year == 2019 & !(scientificName %in% sp2018))
new2021<-subset(d, year == 2021 & !(scientificName %in% sp2019))


sum(d$totalCarabids[d$year == 2017]) #816
sum(d$totalCarabids[d$year == 2016]) #5922
sum(d$totalCarabids[d$year == 2018]) #457
sum(d$totalCarabids[d$year == 2019]) #745
sum(d$totalCarabids[d$year == 2021]) #172
sum(new2017$totalCarabids) #182
sum(new2018$totalCarabids) #177
sum(new2019$totalCarabids) #23
sum(new2021$totalCarabids) #116

length(new2017$scientificName) #17
length(new2018$scientificName) #12
length(new2019$scientificName) #8
length(new2021$scientificName) #5

SRERchange<-data.frame(year = c(2016, 2017, 2018, 2019, 2021), totalCarabids = NA, newSpecies = NA, totalNew = NA)
SRERchange$totalCarabids<-c(5922, 816, 457, 745, 172)
SRERchange$totalNew<-c(NA, 182, 177, 23, 116)
SRERchange$newSpecies<-c(NA, 17, 12, 8, 5)
SRERchange$percentNew<-SRERchange$totalNew / SRERchange$totalCarabids
View(SRERchange)
#repeat this for most and least diverse sites!!



#high diversity site
highdiv<-NeonData$bet_sorting[NeonData$bet_sorting$siteID == 'WOOD',]
highb<-highdiv[highdiv$sampleType == 'carabid',]
highc<-highb %>%
  group_by(year = format(as.Date(setDate),"%Y")) %>%
  summarise(totalCarabids = sum(individualCount), trapsPresent = length(subsampleCode))
View(highc)
highd<- highb %>%
  group_by(scientificName, year = format(as.Date(setDate),"%Y")) %>%
  summarise(totalCarabids = sum(individualCount), trapsPresent = length(subsampleCode))
View(highd)

#high diversity site 2 NOT DONE
highdiv2<-NeonData$bet_sorting[NeonData$bet_sorting$siteID == 'ORNL',]
highb2<-highdiv2[highdiv2$sampleType == 'carabid',]
highc<-highb %>%
  group_by(year = format(as.Date(setDate),"%Y")) %>%
  summarise(totalCarabids = sum(individualCount), trapsPresent = length(subsampleCode))
View(highc)
highd<- highb %>%
  group_by(scientificName, year = format(as.Date(setDate),"%Y")) %>%
  summarise(totalCarabids = sum(individualCount), trapsPresent = length(subsampleCode))
View(highd)

#lots of NAs for this data - what does this mean?
lowdiv<-NeonData$bet_sorting[NeonData$bet_sorting$siteID == 'GUAN',]
lowb<-lowdiv[lowdiv$sampleType == 'carabid',]
lowc<-lowb %>%
  group_by(year = format(as.Date(setDate),"%Y")) %>%
  summarise(totalCarabids = sum(individualCount), trapsPresent = length(subsampleCode))
View(lowc)
lowd<- lowb %>%
  group_by(scientificName, year = format(as.Date(setDate),"%Y")) %>%
  summarise(totalCarabids = sum(individualCount), trapsPresent = length(subsampleCode))
View(lowd)

#another low div site (also 8 species)
lowdiv2<-NeonData$bet_sorting[NeonData$bet_sorting$siteID == 'BARR',]
lowb2<-lowdiv2[lowdiv2$sampleType == 'carabid',]
lowc2<-lowb2 %>%
  group_by(year = format(as.Date(setDate),"%Y")) %>%
  summarise(totalCarabids = sum(individualCount), trapsPresent = length(subsampleCode))
View(lowc2)
lowd2<- lowb2 %>%
  group_by(scientificName, year = format(as.Date(setDate),"%Y")) %>%
  summarise(totalCarabids = sum(individualCount), trapsPresent = length(subsampleCode))
View(lowd2)

#lower div site 3 (17)
lowdiv3<-NeonData$bet_sorting[NeonData$bet_sorting$siteID == 'SJER',]
lowb3<-lowdiv3[lowdiv3$sampleType == 'carabid',]
lowc3<-lowb3 %>%
  group_by(year = format(as.Date(setDate),"%Y")) %>%
  summarise(totalCarabids = sum(individualCount), trapsPresent = length(subsampleCode))
View(lowc3)
lowd3<- lowb3 %>%
  group_by(scientificName, year = format(as.Date(setDate),"%Y")) %>%
  summarise(totalCarabids = sum(individualCount), trapsPresent = length(subsampleCode))
View(lowd3)




write.csv(output,file="C:/Users/avacl/OneDrive/Desktop/speciestotals.csv")







#example for getting total specimens by year
output2<-data.frame(matrix(ncol=10, nrow=0))
colnames(output2)<-c("site",2014,2015,2016,2017,2018,2019,2020,2021,2022)
View(output2)

allsites<-unique(NeonData$bet_expertTaxonomistIDProcessed$siteID)


for (site in allsites){
  print(site)
  sitedata<-NeonData$bet_expertTaxonomistIDProcessed[NeonData$bet_expertTaxonomistIDProcessed$siteID == site,]
  to2014<-length(subset(sitedata, format(as.Date(collectDate), "%Y")==2014)$taxonID)
  to2015<-length(subset(sitedata, format(as.Date(collectDate), "%Y")==2015)$taxonID)
  to2016<-length(subset(sitedata, format(as.Date(collectDate), "%Y")==2016)$taxonID)
  to2017<-length(subset(sitedata, format(as.Date(collectDate), "%Y")==2017)$taxonID)
  to2018<-length(subset(sitedata, format(as.Date(collectDate), "%Y")==2018)$taxonID)
  to2019<-length(subset(sitedata, format(as.Date(collectDate), "%Y")==2019)$taxonID)
  to2020<-length(subset(sitedata, format(as.Date(collectDate), "%Y")==2020)$taxonID)
  to2021<-length(subset(sitedata, format(as.Date(collectDate), "%Y")==2021)$taxonID)
  to2022<-length(subset(sitedata, format(as.Date(collectDate), "%Y")==2022)$taxonID)
  
  newrow<-c(site,to2014,to2015,to2016,to2017,to2018,to2019,to2020,to2021,to2022)
  output2[nrow(output2) +1,]<-newrow
}
View(output2)


