library(lipdR)
library(tidyverse)

#load in list of eDNA lakes:
lakes <- read_csv("lakes.national.era.csv")

#load in LiPD data
load("l380NationalData.Rdata")

#load in new LiPD data
af <- list.files("~/Dropbox/lipdverse/Lakes380National/",full.names = TRUE)

nf <- af[which(file.size(af) < 200000)]

nf <- nf[-which(grepl(nf,pattern = "LakeJohnson_53707.Lakes380"))]

N <- readLipd(nf)

D <- c(D,N)

wl <- data.frame(dataSetName = map_chr(D,\(x) x$dataSetName),
                 Code = map_chr(D,\(x) x$geo$lakeId))

dat <- left_join(lakes,wl)

#missing lakes

missingLakes <- dat |> 
  filter(is.na(dataSetName))

write_csv(missingLakes,"missingLakes.csv")

AD <- D[c(na.omit(dat$dataSetName))]

ts <- extractTs(AD)                 
                 
hsiTs <- ts2tibble(ts) |> 
  filter(paleoData_variableName == "RABD660670")


hsiTs <- rename(hsiTs,Code = geo_lakeId)


#remove old samples
oldDepths <- read_csv("oldsamples.csv") |> 
  group_by(Code) |>
  summarize(depthCut = min(Depth.round)) |> 
  mutate(Code = as.character(Code))

hsiTs <- left_join(hsiTs,oldDepths)
hsiNew <- hsiTs
for(i in 1:nrow(hsiTs)){
  tr <- hsiTs[i,]
  if(!is.na(tr$depthCut)){
    good <- which(tr$depth[[1]] < tr$depthCut)
  }else if(!all(is.na(tr$age[[1]]))){
    good <- which(tr$age[[1]] < 2000)
  }else{
    print("no cut or age")
    good <- seq_along(tr$depth[[1]])
  }
  if(all(is.na(good))){
    print(paste("uh oh",i))
  }
  
  if(!all(is.na(tr$age[[1]]))){
  tr$age[[1]] <- tr$age[[1]][good]
  }
  tr$depth[[1]] <- tr$depth[[1]][good]
  tr$paleoData_values[[1]] <- tr$paleoData_values[[1]][good]
  hsiNew[i,] <- tr
}

toAnalyze <- left_join(lakes, hsiNew) |> 
  filter(!is.na(archiveType))

calcChlIndex <- function(values){
  values[values < 1] <- 1
  vals <- (values - 1) * 1
  return(vals ^ 1)
}

calcEuroIndex <- function(European.Depth,depth){
  return(which(depth <= European.Depth))
}

calcMaoriIndex <- function(European.Depth,Maori.Depth,depth,...){
  return(which(depth <= Maori.Depth & depth > European.Depth))
}

calcPrehumanIndex <- function(Maori.Depth,depth,...){
  return(which(depth > Maori.Depth))
}

calcMeanChlIndex <- function(chl,index){
  return(mean(chl[index],na.rm = TRUE))
}


toAnalyze$maoriIndex = pmap(toAnalyze,calcMaoriIndex)      

#get lake classification
lakeClass <- read_csv("lake.ts.csv") |> 
  mutate(Code = str_remove(str_extract(Lake.Code, "\\.[0-9]+"),"\\.")) |> 
  select(TS,Code)

toAnalyze <- left_join(toAnalyze,lakeClass)

European.Duration <- 2020 - 1830
Maori.Duration <- 1830 - 1260

#calculate windows
toAnalyze <- toAnalyze |> 
  mutate(chlIndex = map(paleoData_values,calcChlIndex),
         euroIndex = map2(European.Depth,depth,calcEuroIndex),
         phIndex = map2(Maori.Depth,depth,calcPrehumanIndex)) |> 
  rowwise() |> 
  mutate(meanEuroChlIndex = calcMeanChlIndex(chlIndex,euroIndex),
         meanMaoriChlIndex = calcMeanChlIndex(chlIndex,maoriIndex),
         meanPhChlIndex = calcMeanChlIndex(chlIndex,phIndex))


#estimate fluxes.
toAnalyze <- toAnalyze |> 
  mutate(meanEuroChlIndexFlux = meanEuroChlIndex * European.Depth / European.Duration , 
         meanMaoriChlIndexFlux = meanMaoriChlIndex * Maori.Depth / Maori.Duration , 
         meanPhChlIndexFlux = meanPhChlIndex * Maori.Depth / Maori.Duration )



#estimate deviations from baseline.
toAnalyze <- toAnalyze |> 
  mutate(euroPct = (meanEuroChlIndex / meanPhChlIndex) * 100 - 100,
         maoriPct = (meanMaoriChlIndex / meanPhChlIndex) * 100 - 100,
         maori2EuroPct = (meanEuroChlIndex / meanMaoriChlIndex) * 100 - 100,
         euroFluxPct = (meanEuroChlIndexFlux / meanPhChlIndexFlux) * 100 - 100,
         maori2EuroFluxPct = (meanEuroChlIndexFlux / meanMaoriChlIndexFlux) * 100 - 100,
         maoriFluxPct = (meanMaoriChlIndexFlux / meanPhChlIndexFlux) * 100 - 100)

#toAnalyze$euroPct[!is.finite(toAnalyze$euroPct)] <- 100

ggplot(toAnalyze) + 
  geom_col(aes( x= Lake,y = euroPct))


ggplot(toAnalyze) + 
  geom_col(aes( x= Lake,y = maoriPct))

#concentration 1st
tal <- pivot_longer(toAnalyze,
                    cols = c("meanEuroChlIndex","meanMaoriChlIndex","meanPhChlIndex"),
                    names_to = "era",
                    names_prefix = "mean",
                    values_to = "meanChlIndex")


tal$era <- str_remove_all(tal$era,"ChlIndex") |> factor(levels = c("Ph","Maori","Euro"))
tal$TS <- factor(tal$TS,levels = c("Oligotrophic","Mesotrophic","Eutrophic","Supertrophic"))

chlIndexMeanTL <- ggplot(tal) + 
  geom_violin(aes(x = era,y = meanChlIndex,fill = era),draw_quantiles = c(.5)) + 
  facet_wrap(TS ~ .,nrow = 1)+
  theme_bw()


ggsave(chlIndexMeanTL,"Chl Index mean by era and trophic level.pdf")


#violin
tal <- pivot_longer(toAnalyze,cols = c("euroPct","maoriPct"),names_to = "era",values_to = "pctChange")
tal$era <- str_remove_all(tal$era,"Pct") |> factor(levels = c("maori","euro"))
tal$TS <- factor(tal$TS,levels = c("Oligotrophic","Mesotrophic","Eutrophic","Supertrophic"))

chlIndexChangeTL <- ggplot(tal) + 
  geom_violin(aes(x = era,y = pctChange,fill = era)) + 
  facet_wrap(TS ~ .,nrow = 1) + 
  ylab("Change in chloropigment concentration (relative to prehuman)") +
  theme_bw()

ggsave(chlIndexChangeTL,"Chl Index change relative to prehuman, by era and trophic level.pdf")



#Maori to Euro change (more lakes)

toAnalyze$TS <- factor(toAnalyze$TS,levels = c("Oligotrophic","Mesotrophic","Eutrophic","Supertrophic"))

chlIndexMaoriToEuroTL <- ggplot(toAnalyze) + 
  geom_violin(aes(x = TS,y = maori2EuroPct,fill = TS)) +
  geom_hline(yintercept = 0,color = "black") +
  ylab("Percent change in chloropigment concentration from Maori to European") +
  theme_bw()

ggsave(chlIndexMaoriToEuroTL,"Chl Index change from Maori to Euro by trophic level.pdf")


flux <- pivot_longer(toAnalyze,cols = c("meanEuroChlIndexFlux","meanMaoriChlIndexFlux"),names_to = "epoch",values_to = "indexFlux")


ggplot(flux) + 
  geom_violin(aes(x = epoch,y = indexFlux))



means <- pivot_longer(toAnalyze,cols = c("meanEuroChlIndex","meanMaoriChlIndex","meanPhChlIndex"),names_to = "epoch",values_to = "meanChl")


violinMeans <- ggplot(means) + 
  geom_violin(aes(x = epoch,y = meanChl)) + 
  scale_x_discrete(labels = c("European","Maori","Pre-human"))+
  ylab("Chloropigments (mean RABD above 1)") + theme_bw()


# maps --------------------------------------------------------------------

library(geoChronR)
res <- toAnalyze
bm <- baseMap(lon = res$geo_longitude, lat = res$geo_latitude,f = .05,projection = "mercator") +
  theme_bw() +
  scale_x_continuous(name = NULL,expand = c(0,0)) +
  scale_y_continuous(name = NULL,expand = c(0,0))

climits <- c(-60,600)


euroConcChange <- bm +  geom_point(data = toAnalyze,aes(x = geo_longitude, 
                                                        y = geo_latitude), color = "black"
                                   , size = 5) +
  geom_point(data = toAnalyze,aes(x = geo_longitude, 
                                  y = geo_latitude, 
                                  color = euroPct) , size = 4) +
  scale_color_gradient2("Percent change",limits = climits,low = "#543005",mid = "white",high = "DarkGreen",transform = scales::pseudo_log_trans(sigma = 3,base = exp(1)),breaks = c(-40,-10,0,10,40,100,400)) + 
  ggtitle("European era chloropigment concentration\n (relative to prehuman)") +
  theme(legend.position = "inside",legend.position.inside = c(.8,.2))



maoriConcChange <- bm +  geom_point(data = toAnalyze,aes(x = geo_longitude, 
                                                         y = geo_latitude), color = "black"
                                    , size = 5) +
  geom_point(data = toAnalyze,aes(x = geo_longitude, 
                                  y = geo_latitude, 
                                  color = maoriPct) , size = 4) +
  scale_color_gradient2("Percent change",limits = climits,low = "#543005",mid = "white",high = "DarkGreen",transform = scales::pseudo_log_trans(sigma = 3,base = exp(1)),breaks = c(-40,-10,0,10,40,100,400)) + 
  ggtitle("Maori era chloropigment concentration\n (relative to prehuman)") +
  theme(legend.position = "inside",legend.position.inside = c(.8,.2))




deltaConcChange <- bm +  geom_point(data = toAnalyze,aes(x = geo_longitude, 
                                                         y = geo_latitude), color = "black"
                                    , size = 5) +
  geom_point(data = toAnalyze,aes(x = geo_longitude, 
                                  y = geo_latitude, 
                                  color = euroPct - maoriPct) , size = 4) +
  scale_color_gradient2("Delta percent",low = "#543005",mid = "white",high = "DarkGreen",transform = scales::pseudo_log_trans(sigma = 2,base = exp(1)), breaks = c( -30, -10, 0, 10, 30, 100, 300)) + 
  ggtitle("Difference between European and Maori impacts on Concentration") +
  theme(legend.position = "inside",legend.position.inside = c(.8,.2))


library(egg)
out <- ggarrange(plots = list(euroConcChange,maoriConcChange,deltaConcChange),nrow = 1,)

ggsave(plot = out, filename = "percent change maps.pdf")



# flux pct change maps ----------------------------------------------------



climits <- c(range(c(toAnalyze$euroFluxPct,toAnalyze$maoriFluxPct),na.rm = TRUE))


euroFluxChange <- bm +  geom_point(data = toAnalyze,aes(x = geo_longitude, 
                                                        y = geo_latitude), color = "black"
                                   , size = 5) +
  geom_point(data = toAnalyze,aes(x = geo_longitude, 
                                  y = geo_latitude, 
                                  color = euroFluxPct) , size = 4) +
  scale_color_gradient2("Percent change",limits = climits, low = "#543005",mid = "white",high = "DarkGreen",transform = scales::pseudo_log_trans(sigma = 3,base = exp(1)),breaks = c(-40,-10,0,10,40,100,400)) + 
  ggtitle("European era chloropigment flux\n (relative to prehuman)") +
  theme(legend.position = "inside",legend.position.inside = c(.8,.2))



maoriFluxChange <- bm +  geom_point(data = toAnalyze,aes(x = geo_longitude, 
                                                         y = geo_latitude), color = "black"
                                    , size = 5) +
  geom_point(data = toAnalyze,aes(x = geo_longitude, 
                                  y = geo_latitude, 
                                  color = maoriFluxPct) , size = 4) +
  scale_color_gradient2("Percent change",limits = climits,low = "#543005",mid = "white",high = "DarkGreen",transform = scales::pseudo_log_trans(sigma = 3,base = exp(1)),breaks = c(-40,-10,0,10,40,100,400)) + 
  ggtitle("Maori era chloropigment flux\n (relative to prehuman)") +
  theme(legend.position = "inside",legend.position.inside = c(.8,.2))




deltaFluxChange <- bm +  geom_point(data = toAnalyze,aes(x = geo_longitude, 
                                                         y = geo_latitude), color = "black"
                                    , size = 5) +
  geom_point(data = toAnalyze,aes(x = geo_longitude, 
                                  y = geo_latitude, 
                                  color = euroFluxPct - maoriFluxPct) , size = 4) +
  scale_color_gradient2("Delta percent",low = "#543005",mid = "white",high = "DarkGreen",transform = scales::pseudo_log_trans(sigma = 2,base = exp(1)), breaks = c( -30, -10, 0, 10, 30, 100, 300)) + 
  ggtitle("Difference between European and\n Maori impacts on flux") +
  theme(legend.position = "inside",legend.position.inside = c(.8,.2))


library(egg)
out <- ggarrange(plots = list(euroFluxChange,maoriFluxChange,deltaFluxChange),nrow = 1,)

ggsave(plot = out, filename = "flux percent change maps.pdf",width = 13)


