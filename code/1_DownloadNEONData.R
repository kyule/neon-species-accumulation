### This is the  file for downloading NEON data

#### BEFORE RUNNING

#### Users must define their own paths and NEON token

# datapath<-"user defined path"
# codepath<-"user defined path"
# Neon_Token<-"user's token"

# use the github version for now
#remotes::install_github("NEONScience/NEON-utilities/neonUtilities")

#### And users should define their variables for the NEON data download

# data product of interest, below is for Carabids
product<-"DP1.10022.001"

# start and end dates, must be formatted as "YYYY-MM" or NA for all time
start<-"2012-01"
end<-"2024-05"

# sites, concatenated list of siteCodes or use "all" for all sites
sites<-"all"

# And user should remove below
source("/Users/kelsey/Github/neon-species-accumulation/configini.R")

# User can change other options in the loadByProduct function as desired

# load package
library(neonUtilities)


# Use loadByProduct to download the data
NeonData<-loadByProduct(dpID=product,
                          site=sites,
                          startdate=start,
                          enddate=end,
                          check.size=FALSE, 
                          release= 'RELEASE-2025',
                          include.provisional=FALSE,
                          token=Neon_Token)
save(NeonData, file=paste0(datapath,"NeonData.Robj"))
