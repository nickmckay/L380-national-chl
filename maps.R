source("calibrationFunctions.R")

library(lipdR)
library(geoChronR)
library(tidyverse)
D <- readLipd("~/Dropbox/lipdverse/Lakes380National/")

ts <- extractTs(D) %>% ts2tibble()

culturalEpoch <- function(year){
  culturalEpoch <- rep(NA,times = length(year))
  culturalEpoch[which(year >= 1840)] <- "European"
  culturalEpoch[which(year < 1840 & year >= 1340)] <- "Maori"
  culturalEpoch[which(year < 1340)] <- "Pre-human"
  return(culturalEpoch)
}

chlgood <- filter(ts,paleoData_variableName == "RABD660670" & paleoData_useInL380National == "TRUE")

#calculate means and differences. 

getMeansAndDeltas <- function(sumStats,...){
  
  phi <- which(sumStats$culturalEpoch == "Pre-human")
  mi <- which(sumStats$culturalEpoch == "Maori")
  ei <- which(sumStats$culturalEpoch == "European")
  
  
  out <- tibble(prehumanConcMean = sumStats$concMean[phi],
                maoriConcMean = sumStats$concMean[mi], 
                euroConcMean = sumStats$concMean[ei], 
                prehumanFluxMean = sumStats$fluxMean[phi],
                maoriFluxMean = sumStats$fluxMean[mi], 
                euroFluxMean = sumStats$fluxMean[ei], 
                maoriConcDelta = sumStats$concDelta[phi], 
                euroConcDelta = sumStats$concDelta[mi], 
                maoriFluxDelta = sumStats$fluxDelta[phi], 
                euroFluxDelta = sumStats$fluxDelta[mi], 
                maoriFluxSdDelta = sumStats$fluxDeltaAsSd[phi], 
                euroFluxSdDelta = sumStats$fluxDeltaAsSd[mi], 
                maoriFluxTTDelta = sumStats$flux.ttest.p[phi], 
                euroFluxTTDelta = sumStats$flux.ttest.p[mi] 

               )
  
  return(out)
}


chlgood$sumStats <- vector(mode = "list",length = nrow(chlgood))
for(i in 1:nrow(chlgood)){
  year <- convertBP2AD(chlgood$age[[i]])
  thisLake <- tibble(year = year,
                   depth = chlgood$depth[[i]],
                   RABD660670 = chlgood$paleoData_values[[i]]) %>% 
    filter(!duplicated(year), year > 500) %>% 
    mutate(culturalEpoch = culturalEpoch(year),
           dbd = estBulkDensity(depth),
           deltaYear = abs(deltas(year)),
           deltaDepth = abs(deltas(depth)),
           sedRate = deltaDepth/deltaYear,
           chlEst = calcConc(RABD660670),
           chlFlux = chlEst * dbd * sedRate) %>% 
    filter(depth > 1,
           chlFlux < 1e6) 
  
  
  sumStats <- thisLake %>% 
    group_by(culturalEpoch) %>% 
    summarize(concMean = mean(chlEst,na.rm = TRUE),
              concSd = sd(chlEst,na.rm = TRUE),
              fluxMean = mean(chlFlux,na.rm = TRUE),
              fluxSd = sd(chlFlux,na.rm = TRUE)) %>% 
    mutate(concDelta = c(NA,-diff(concMean)),
           fluxDelta = c(NA,-diff(fluxMean)),
           concDeltaAsSd = concDelta/concSd,
           fluxDeltaAsSd = fluxDelta/fluxSd)
  
  #t.test
  ec <- thisLake$chlEst[thisLake$culturalEpoch == "European"]
  mc <- thisLake$chlEst[thisLake$culturalEpoch == "Maori"]
  pc <- thisLake$chlEst[thisLake$culturalEpoch == "Pre-human"]
  ef <- thisLake$chlFlux[thisLake$culturalEpoch == "European"]
  mf <- thisLake$chlFlux[thisLake$culturalEpoch == "Maori"]
  pf <- thisLake$chlFlux[thisLake$culturalEpoch == "Pre-human"]
  
  emcd <- t.test(ec,mc)$p.value
  mpcd <- t.test(pc,mc)$p.value
  emfd <- t.test(ef,mf)$p.value
  mpfd <- t.test(pf,mf)$p.value
  
  sumStats$conc.ttest.p <- c(NA, emcd,mpcd)
  sumStats$flux.ttest.p <- c(NA, emfd,mpfd)
  
  chlgood$sumStats[[i]] <- sumStats
  
}



res <- pmap_dfr(chlgood,getMeansAndDeltas) %>% 
  bind_cols(chlgood) %>% 
  select(starts_with("geo"),ends_with("Delta"),ends_with("Mean")) %>% 
  mutate(across(.cols = everything(),.fns = ~ na_if(.x,"NaN"))) %>% 
  mutate(maoriFluxSignificant = res$maoriFluxTTDelta < 0.01,
         euroFluxSignificant = res$euroFluxTTDelta < 0.01)



dplotFlux <- ggplot() + 
  geom_density(aes(x = res$maoriFluxDelta, fill = "Maori arrival"), alpha = 0.5) + 
  geom_density(aes(x = res$euroFluxDelta, fill = "European arrival"), alpha = 0.5) +
  geom_vline(aes(xintercept = 0)) +
  xlab("Change in mean Chla flux betwen epochs")+
  ggtitle("Chla flux shift with epoch") +
  scale_fill_brewer("Cultural Epoch",palette = "Dark2") +
  theme_bw()

dplotConc <- ggplot() + 
  geom_density(aes(x = res$maoriConcDelta, fill = "Maori arrival"), alpha = 0.5) + 
  geom_density(aes(x = res$euroConcDelta, fill = "European arrival"), alpha = 0.5) +
  geom_vline(aes(xintercept = 0)) +
  xlab("Change in mean Chla concentration between epochs")+
  ggtitle("Chla concentration shift with epoch") +
  scale_fill_brewer("Cultural Epoch",palette = "Dark2") +
  theme_bw()

library(geoChronR)
bm <- baseMap(lon = res$geo_longitude, lat = res$geo_latitude,f = .05,projection = "mercator") +
  theme_bw() +
  scale_x_continuous(name = NULL,expand = c(0,0)) +
  scale_y_continuous(name = NULL,expand = c(0,0))

climits <- c(range(c(res$maoriConcDelta,res$euroConcDelta),na.rm = TRUE))

td1 <- drop_na(res,maoriConcDelta) %>% 
  arrange(abs(maoriConcDelta))


MaoriMapOwnScale <- bm +  geom_point(data = td1,aes(x = geo_longitude, 
                                            y = geo_latitude), color = "black"
                             , size = 6) +
  geom_point(data = td1,aes(x = geo_longitude, 
                            y = geo_latitude, 
                            color = maoriConcDelta) , size = 5) +
  scale_color_gradient2("Change in Chla concentration",low = "#543005",mid = "white",high = "Dark Green") + 
  ggtitle("Impact of Maori Arrival on Concentration")

MaoriMap <- bm +  geom_point(data = td1,aes(x = geo_longitude, 
                                            y = geo_latitude), color = "black"
                             , size = 6) +
  geom_point(data = td1,aes(x = geo_longitude, 
                            y = geo_latitude, 
                            color = maoriConcDelta) , size = 5) +
  scale_color_gradient2("Change in Chla concentration",low = "#543005",mid = "white",high = "Dark Green",limits = climits) + 
  ggtitle("Impact of Maori Arrival on Concentration")


td2 <- drop_na(res,euroConcDelta) %>% 
  arrange(abs(euroConcDelta))

EuroMap <- bm +  geom_point(data = td2,aes(x = geo_longitude, 
                                           y = geo_latitude), color = "black"
                            , size = 6) +
  geom_point(data = td2,aes(x = geo_longitude, 
                            y = geo_latitude, 
                            color = euroConcDelta) , size = 5) +
  scale_color_gradient2("Change in Chla concentration",low = "#543005",mid = "white",high = "Dark Green",limits = climits) + 
  ggtitle("Impact of European Arrival on Chl")

allMaps <- list(dplotConc,MaoriMapOwnScale,MaoriMap,EuroMap)


#gridExtra::grid.arrange(grobs = allMaps,nrow = 2,... = )
mapOut <- egg::ggarrange(plots = allMaps,nrow = 2)
ggsave(filename = "RABD shift conc - first look.pdf",mapOut,height = 15, width = 15)




# repeat with flux --------------------------------------------------------
climits <- c(range(c(res$maoriFluxDelta,res$euroFluxDelta),na.rm = TRUE))

td1 <- drop_na(res,maoriFluxDelta) %>% 
  arrange(abs(maoriFluxDelta))


MaoriMapOwnScale <- bm +  geom_point(data = td1,aes(x = geo_longitude, 
                                                    y = geo_latitude), color = "black"
                                     , size = 6) +
  geom_point(data = td1,aes(x = geo_longitude, 
                            y = geo_latitude, 
                            color = maoriFluxDelta) , size = 5) +
  scale_color_gradient2("Change in Chla Flux",low = "#543005",mid = "white",high = "Dark Green") + 
  ggtitle("Impact of Maori Arrival on Flux")

MaoriMap <- bm +  geom_point(data = td1,aes(x = geo_longitude, 
                                            y = geo_latitude), color = "black"
                             , size = 6) +
  geom_point(data = td1,aes(x = geo_longitude, 
                            y = geo_latitude, 
                            color = maoriFluxDelta) , size = 5) +
  scale_color_gradient2("Change in Chla Flux",low = "#543005",mid = "white",high = "Dark Green",limits = climits) + 
  ggtitle("Impact of Maori Arrival on Flux")


td2 <- drop_na(res,euroFluxDelta) %>% 
  arrange(abs(euroFluxDelta))

EuroMap <- bm +  geom_point(data = td2,aes(x = geo_longitude, 
                                           y = geo_latitude), color = "black"
                            , size = 6) +
  geom_point(data = td2,aes(x = geo_longitude, 
                            y = geo_latitude, 
                            color = euroFluxDelta) , size = 5) +
  scale_color_gradient2("Change in Chla Flux",low = "#543005",mid = "white",high = "Dark Green",limits = climits) + 
  ggtitle("Impact of European Arrival on Chl Flux")

allMaps <- list(dplotFlux,MaoriMapOwnScale,MaoriMap,EuroMap)


#gridExtra::grid.arrange(grobs = allMaps,nrow = 2,... = )
mapOut <- egg::ggarrange(plots = allMaps,nrow = 2)
ggsave(filename = "RABD shift flux - first look.pdf",mapOut,height = 15, width = 15)



#plot vs elevation?
# 
# ggplot(res) + geom_point(aes(x = modernDel, y = geo_elevation))
# 
# ggplot(res) + geom_point(aes(x = euroDelta, y = geo_elevation))
# ggplot(res) + geom_point(aes(x = maoriDelta, y = geo_elevation))
# 
# 
# ggplot(resFlux) + geom_point(aes(x = modernDel, y = geo_elevation))
# 
# ggplot(resFlux) + geom_point(aes(x = euroDelta, y = geo_elevation))
# ggplot(resFlux) + geom_point(aes(x = maoriDelta, y = geo_elevation))


# Concentration amount maps -----------------------------------------------

legTitle <- expression(paste("Chloropigment concentration (",mu,"g ", cm^-3,")"))

climits <- range(c(res$prehumanConcMean, res$maoriConcMean,res$euroConcMean),na.rm = TRUE)

td0 <- drop_na(res,prehumanConcMean) %>% 
  arrange(abs(prehumanConcMean))

PreHuman <- bm +  geom_point(data = td0,aes(x = geo_longitude, 
                                            y = geo_latitude), color = "black"
                             , size = 6) +
  geom_point(data = td0,aes(x = geo_longitude, 
                            y = geo_latitude, 
                            color = prehumanConcMean) , size = 5) +
  scale_color_viridis_c(legTitle,limits = climits) + 
  ggtitle("Mean chla pre-Maori arrival")

PreHumanOwnScale <- bm +  geom_point(data = td0,aes(x = geo_longitude, 
                                            y = geo_latitude), color = "black"
                             , size = 6) +
  geom_point(data = td0,aes(x = geo_longitude, 
                            y = geo_latitude, 
                            color = prehumanConcMean) , size = 5) +
  scale_color_viridis_c(legTitle) + 
  ggtitle("Mean chla pre-Maori arrival")


td1 <- drop_na(res,maoriConcMean) %>% 
  arrange(abs(maoriConcMean))

MaoriMap <- bm +  geom_point(data = td1,aes(x = geo_longitude, 
                                            y = geo_latitude), color = "black"
                             , size = 6) +
  geom_point(data = td1,aes(x = geo_longitude, 
                            y = geo_latitude, 
                            color = maoriConcMean) , size = 5) +
  scale_color_viridis_c(legTitle,limits = climits) + 
  ggtitle("Maori")


td2 <- drop_na(res,euroConcMean) %>% 
  arrange(abs(euroConcMean))

EuroMap <- bm +  geom_point(data = td2,aes(x = geo_longitude, 
                                           y = geo_latitude), color = "black"
                            , size = 6) +
  geom_point(data = td2,aes(x = geo_longitude, 
                            y = geo_latitude, 
                            color = euroConcMean) , size = 5) +
  scale_color_viridis_c(legTitle,limits = climits) + 
  ggtitle("European")


allMaps <- list(PreHumanOwnScale,PreHuman,MaoriMap,EuroMap)


#gridExtra::grid.arrange(grobs = allMaps,nrow = 2,... = )
mapOut <- egg::ggarrange(plots = allMaps,nrow = 2)
ggsave(filename = "Chla Concentration Maps.pdf",mapOut,height = 15, width = 15)


ggplot(res)+geom_point(aes(x = prehumanConcMean,y = geo_elevation))+geom_point(aes(x = euroConcMean,y = geo_elevation),color = "red")


# Flux amount maps -----------------------------------------------
#constant scale

climits <- range(c(res$prehumanFluxMean, res$maoriFluxMean,res$euroFluxMean),na.rm = TRUE)

climits[1] <- 0

td0 <- drop_na(res,prehumanFluxMean) %>% 
  arrange(abs(prehumanFluxMean))

legTitle <- expression(paste("Chloropigment flux (",mu,"g ", cm^-2," ",yr^-1,")"))

PreHuman <- bm +  geom_point(data = td0,aes(x = geo_longitude, 
                                            y = geo_latitude), color = "black"
                             , size = 6) +
  geom_point(data = td0,aes(x = geo_longitude, 
                            y = geo_latitude, 
                            color = prehumanFluxMean) , size = 5) +
  scale_color_viridis_c(limits = climits) + 
  labs(color = legTitle) +
  ggtitle("Pre-human arrival")

PreHumanOwnScale <- bm +  geom_point(data = td0,aes(x = geo_longitude, 
                                            y = geo_latitude), color = "black"
                             , size = 6) +
  geom_point(data = td0,aes(x = geo_longitude, 
                            y = geo_latitude, 
                            color = prehumanFluxMean) , size = 5) +
  scale_color_viridis_c() + 
  labs(color = legTitle) +
  ggtitle("Pre-human arrival (separate scale)")


td1 <- drop_na(res,maoriFluxMean) %>% 
  arrange(abs(maoriFluxMean))

MaoriMap <- bm +  geom_point(data = td1,aes(x = geo_longitude, 
                                            y = geo_latitude), color = "black"
                             , size = 6) +
  geom_point(data = td1,aes(x = geo_longitude, 
                            y = geo_latitude, 
                            color = maoriFluxMean) , size = 5) +
  scale_color_viridis_c(limits = climits) + 
  labs(color = legTitle) +
  
  ggtitle("Maori")


td2 <- drop_na(res,euroFluxMean) %>% 
  arrange(abs(euroFluxMean))

EuroMap <- bm +  geom_point(data = td2,aes(x = geo_longitude, 
                                           y = geo_latitude), color = "black"
                            , size = 6) +
  geom_point(data = td2,aes(x = geo_longitude, 
                            y = geo_latitude, 
                            color = euroFluxMean) , size = 5) +
  scale_color_viridis_c(limits = climits) + 
  labs(color = legTitle) +
  
  ggtitle("European")


allMaps <- list(PreHumanOwnScale,PreHuman,MaoriMap,EuroMap)

#allMaps <- list(PreHuman,MaoriMap,EuroMap, ModernMap)


#gridExtra::grid.arrange(grobs = allMaps,nrow = 2,... = )
mapOut <- egg::ggarrange(plots = allMaps,nrow = 2)
ggsave(filename = "Chla flux Maps - constant scale.pdf",mapOut,height = 15, width = 15)



# Flux maps changing scale ------------------------------------------------




