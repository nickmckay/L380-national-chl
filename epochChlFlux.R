#epoch flux calculations

source("calibrationFunctions.R")

library(lipdR)
library(geoChronR)
library(tidyverse)

D <- readLipd("~/Dropbox/lipdverse/Lakes380National/")

ts <- extractTs(D) %>% ts2tibble()

chlgood <- filter(ts,paleoData_variableName == "RABD660670" & paleoData_useInL380National == "TRUE")

legTitle <- expression(paste("Chloropigment flux (",mu,"g ", cm^-2," ",yr^-1,")"))
concTitle <- expression(paste("Chloropigment concentration (",mu,"g ", g^-1,")"))

for(i in 1:nrow(chlgood)){
  year <- convertBP2AD(chlgood$age[[i]])
  toplot <- tibble(year = year,
                   depth = chlgood$depth[[i]],
                   RABD660670 = chlgood$paleoData_values[[i]]) %>% 
    filter(!duplicated(year)) %>% 
    mutate(culturalEpoch = culturalEpoch(year),
           dbd = estBulkDensity(depth),
           deltaYear = abs(deltas(year)),
           deltaDepth = abs(deltas(depth)),
           sedRate = deltaDepth/deltaYear,
           chlEst = calcConc(RABD660670),
           chlFlux = chlEst * dbd * sedRate) %>% 
    filter(depth > 1,
           chlFlux < 1e6)
  
  fo <- ggplot(toplot) + 
    geom_line(aes(x = year, y = chlFlux,color = culturalEpoch)) +
    scale_color_brewer("Cultural Epoch",palette = "Dark2") + 
    labs(x = "Year AD",y = "Chloropigment flux (ug/cm2/yr)",title = chlgood$geo_siteName[[i]]) +
    theme_bw()+
    theme(legend.position = c(.2,.8))
  
  ggsave(plot = fo,filename = glue::glue("plots/ChlFlux/{chlgood$geo_siteName[[i]]}-ChlFlux-CulturalEpoch.pdf"),width = 6,height = 4)
  
  fo2 <- fo + xlim(c(1,2020))
  ggsave(plot = fo2,filename = glue::glue("plots/ChlFlux-2k/{chlgood$geo_siteName[[i]]}-ChlFlux-CulturalEpoch.pdf"),width = 6,height = 4)
  
  fo <- ggplot(toplot) + 
    geom_line(aes(x = year, y = chlEst,color = culturalEpoch)) +
    scale_color_brewer("Cultural Epoch",palette = "Dark2") + 
    labs(x = "Year AD",y = "Chloropigment concentration (ug/cm3)",title = chlgood$geo_siteName[[i]]) +
    theme_bw()+
    theme(legend.position = c(.2,.8))
  
  ggsave(plot = fo,filename = glue::glue("plots/ChlConc/{chlgood$geo_siteName[[i]]}-ChlFlux-CulturalEpoch.pdf"),width = 6,height = 4)
  
  

# Violin plots ------------------------------------------------------------

  concViolin <- ggplot(toplot) + 
      geom_violin(aes(x = culturalEpoch, y = chlEst, fill = culturalEpoch)) +
      ggtitle(glue::glue("{chlgood$geo_siteName[[i]]} - Concentration")) +
      scale_fill_brewer("Cultural Epoch",palette = "Dark2") +
      theme_bw()+
      ylab(concTitle) +
      xlab("")  
  
  ggsave(plot = concViolin,filename = glue::glue("plots/ChlConcViolin/{chlgood$geo_siteName[[i]]}.pdf"),width = 6,height = 4)
  
  fluxViolin <- ggplot(toplot) + 
    geom_violin(aes(x = culturalEpoch, y = chlFlux, fill = culturalEpoch)) +
    ggtitle(glue::glue("{chlgood$geo_siteName[[i]]} - Flux")) +
    scale_fill_brewer("Cultural Epoch",palette = "Dark2") +
    theme_bw()+
    ylab(legTitle) +
    xlab("")  
  
  ggsave(plot = fluxViolin,filename = glue::glue("plots/ChlFluxViolin/{chlgood$geo_siteName[[i]]}.pdf"),width = 6,height = 4)
  
  
}

ad <- list.dirs("plots")
for(d in 2:length(ad)){
  tp <- file.path(getwd(),ad[d])
system(glue::glue("rm {file.path(tp,'allplots.pdf')}"))
system(glue::glue("cd {file.path(tp)};pdftk *.pdf cat output allplots.pdf"))
}
