
data(spider)
spider.next<-iNEXT(spider,q=c(0,1,2),datatype = "abundance")

spider.next$DataInfo
spider.next$iNextEst
spider.next$AsyEst

data(ant)
ant.next<-iNEXT(ant,q=c(0,1,2),datatype = "incidence_freq")

ant.next$DataInfo
ant.next$iNextEst
ant.next$AsyEst



estimateD(ant, datatype="incidence_freq", base="coverage", level= 0.985, conf=0.95)

data(ciliates)
cil.next<-iNEXT(ciliates,q=c(0,1,2),datatype = "incidence_raw")

cil.next$DataInfo
cil.next$iNextEst
cil.next$AsyEst













# Stable version
install.packages('rarestR')
# Development version
remotes::install_github('pzhaonet/rarestR')
library(rarestR)
data("share")
es(share, m = 80)
es(share, method = "b", m = 100)
ess(share)
ess(share, m = 100)
ess(share, m = 100, index = "ESS")
Output_tes <- tes(share[1,])
Output_tes
plot(Output_tes)
Output_tess <- tess(share[1:2,])
Output_tess
plot(Output_tess)

data(bird)
## Not run:
out <- iNEXT(bird, datatype="abundance")
ggiNEXT(out)

ChaoRichness(bird)

out <- iNEXT(ciliates, datatype = "incidence_raw")
ggiNEXT(out)


test<-comm %>% select(-year)
es(test,m=100)
es(test, method = "b", m = 100)
ess(test)

Output_tess <- tess(test[1:2,])
Output_tess

plot(Output_tess)

library(iNEXT)
test<-test[,-41]

test2<-t(test)
test2[which(test2>0)]<-1

test3<-t(comm[,-19])
test3[which(test3>0)]<-1

test4<-t(comm[,-21])
test4[which(test4>0)]<-1

testa<-list(bart=test2,guan=test3,puum=test4)


#a<-iNEXT(test, q=0, datatype='abundance', size=NULL, endpoint=NULL, knots=40, se=TRUE, conf=0.95, nboot=50)
b<-iNEXT(testa,q=0,datatype='incidence_raw')
ggiNEXT(b, type=3, facet.var="None")
ggiNEXT(b, type=2, facet.var="None")
ggiNEXT(b, type=1, facet.var="None")

b$AsyEst->c

#estimateD(testa)

lapply(ciliates, as.incfreq)->a


#iNEXT: https://github.com/AnneChao/iNEXT