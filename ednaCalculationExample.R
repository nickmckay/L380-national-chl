#calculate shannon diversity index

library(vegan)
library(purrr)
library(lipdR)
library(tidyverse)
library(geoChronR) #remotes::install_github("nickmckay/geoChronR")

#D <- readLipd("path/to/national/data")
#or
load("l380NationalData.Rdata")


#function to extract table from LiPD and calculate diversity index
calcDiv <- function(L){
  #find eDNA data
  tabs <- getMeasurementTables(L,pc = "paleo")  
  
  wt <- which.max(map_dbl(tabs,ncol))
  
  #silly way to get the paleo and table for later use
  nums <- str_extract_all(names(wt),pattern = "\\d{1,}")[[1]]
  p <- as.numeric(nums[1])
  m <- as.numeric(nums[2])
  
  #check for proxy in table
  allProxy <- map(L$paleoData[[p]]$measurementTable[[m]],safely(pluck,otherwise = "none"),"proxy") %>% 
    unlist() %>% unique()
  
  if(!"eDNA" %in% allProxy){#if there are no eDNA data, return an unmodified lipd object
    return(L)
  }#otherwise, keep going
  
  #isolate the data and get a diversity index
  ednaTable <- tabs[[which.max(map_dbl(tabs,ncol))]]
  asvsCols <- which(map_dbl(names(ednaTable),str_length) > 100)
  asvOnly <- select(ednaTable,asvsCols) %>% 
    transmute(across(everything(),as.numeric)) %>% 
    as.matrix()
  div <- vegan::diversity(asvOnly,index = "shannon")
  
  L <- lipdR::createColumn(L,
               paleo.or.chron.number = p,
               table.number = m,
               variableName = "diversityIndex",
               units = "unitless",
               values = div,
               additional.metadata = list(method = "shannon index calculated by vegan")
                )
  
  return(L)
}

#calculate diversity for for all
Dd <- map(D,calcDiv)

#
tsd <- extractTs(Dd) %>% ts2tibble()

#grab only diversity data
tsDiv <- dplyr::filter(tsd,paleoData_variableName == "diversityIndex")

#fucntion to assign cultural epoch
culturalEpoch <- function(year){
  culturalEpoch <- rep(NA,times = length(year))
  culturalEpoch[which(year >= 1840)] <- "European"
  culturalEpoch[which(year < 1840 & year >= 1340)] <- "Maori"
  culturalEpoch[which(year < 1340)] <- "Pre-human"
  return(culturalEpoch)
}


#calculate summary statistics by epoch
tsDiv$sumStats <- vector(mode = "list",length = nrow(chlgood))
for(i in 1:nrow(tsDiv)){
  year <- convertBP2AD(tsDiv$age[[i]])
  thisLake <- tibble(year = year,
                     depth = tsDiv$depth[[i]],
                     diversity = tsDiv$paleoData_values[[i]]) %>% 
    filter(!duplicated(year), year > 500) %>% 
    mutate(culturalEpoch = culturalEpoch(year)) %>% 
    filter(depth > 0) 
  
  
  sumStats <- thisLake %>% 
    group_by(culturalEpoch) %>% 
    summarize(divMean = mean(diversity,na.rm = TRUE),
              divSd = sd(diversity,na.rm = TRUE)) %>% 
    mutate(divDelta = c(NA,-diff(divMean)),
           divDeltaAsSd = divDelta/divSd)
  

  tsDiv$sumStats[[i]] <- sumStats
  
}


#function to flatten results into a row in the tibble
getMeansAndDeltas <- function(sumStats,...){
  
  phi <- which(sumStats$culturalEpoch == "Pre-human")
  mi <- which(sumStats$culturalEpoch == "Maori")
  ei <- which(sumStats$culturalEpoch == "European")
  
  
  out <- tibble(prehumanDivMean = sumStats$divMean[phi],
                maoriDivMean = sumStats$divMean[mi], 
                euroDivMean = sumStats$divMean[ei], 
                maoriDivDelta = sumStats$divDelta[phi], 
                euroDivDelta = sumStats$divDelta[mi], 
                maoriDivSdDelta = sumStats$divDeltaAsSd[phi], 
                euroDivSdDelta = sumStats$divDeltaAsSd[mi]
  )
  
  return(out)
}

#only take those with data in pre human, maori and euro epochs
alll <- map_dbl(tsDiv$sumStats,~ nrow(.x))
tsDivGood <- filter(tsDiv,alll == 3)

#create table for plotting
res <- pmap_dfr(tsDivGood,getMeansAndDeltas) %>% 
  bind_cols(tsDivGood) %>% 
  select(starts_with("geo"),ends_with("Delta"),ends_with("Mean")) %>% 
  mutate(across(.cols = everything(),.fns = ~ na_if(.x,"NaN")))



#create the basemap

bm <- geoChronR::baseMap(lon = res$geo_longitude, lat = res$geo_latitude,f = .05,projection = "mercator") +
  theme_bw() +
  scale_x_continuous(name = NULL,expand = c(0,0)) +
  scale_y_continuous(name = NULL,expand = c(0,0))


climits <- c(range(c(res$prehumanDivMean,res$euroDivMean),na.rm = TRUE))


#pre human mean
phMean <- bm +  geom_point(data = td1,aes(x = geo_longitude, 
                                y = geo_latitude), color = "black"
                 , size = 6) +
  geom_point(data = td1,aes(x = geo_longitude, 
                            y = geo_latitude, 
                            color = prehumanDivMean) , size = 5) +
  scale_color_distiller(guide = "colorbar","eDNA diversity",direction = 1, palette = "Purples",limits = climits) + 
  ggtitle("Pre-human mean eDNA diversity")

#euro mean
euroMean <- bm +  geom_point(data = td1,aes(x = geo_longitude, 
                                y = geo_latitude), color = "black"
                 , size = 6) +
  geom_point(data = td1,aes(x = geo_longitude, 
                            y = geo_latitude, 
                            color = euroDivMean) , size = 5) +
  scale_color_distiller(guide = "colorbar","eDNA diversity",direction = 1, palette = "Purples",limits = climits) + 
  ggtitle("European epoch mean eDNA diversity")





#map delta units

climits <- c(range(c(res$maoriDivDelta,res$euroDivDelta),na.rm = TRUE))
climits <- c(-max(abs(climits)),max(abs(climits)))

td1 <- drop_na(res,maoriDivDelta) %>% 
  arrange(abs(maoriDivDelta))


MaoriMap <- bm +  geom_point(data = td1,aes(x = geo_longitude, 
                                                    y = geo_latitude), color = "black"
                                     , size = 6) +
  geom_point(data = td1,aes(x = geo_longitude, 
                            y = geo_latitude, 
                            color = maoriDivDelta) , size = 5) +
  scale_color_distiller(guide = "colorbar","Change in eDNA diversity",palette = "PuOr",limits = climits,direction = 1) + 
  ggtitle("Impact of Maori Arrival on eDNA diversity")



EuroMap <- bm +  geom_point(data = td1,aes(x = geo_longitude, 
                                                    y = geo_latitude), color = "black"
                                     , size = 6) +
  geom_point(data = td1,aes(x = geo_longitude, 
                            y = geo_latitude, 
                            color = euroDivDelta) , size = 5) +
  scale_color_distiller(guide = "colorbar","Change in eDNA diversity",palette = "PuOr",limits = climits,direction = 1) + 
  ggtitle("Impact of European Arrival on eDNA diversity")


#plot the maps together.

egg::ggarrange(plots = list(phMean, euroMean,MaoriMap,EuroMap),ncol = 2)

