#quick epoch plots

library(lipdR)
library(geoChronR)
library(tidyverse)

D <- readLipd("~/Dropbox/lipdverse/Lakes380National/")

ts <- extractTs(D) %>% ts2tibble()

culturalEpoch <- function(year){
  culturalEpoch <- matrix(NA,nrow = length(year))
  culturalEpoch[which(year >= 1840)] <- "European"
  culturalEpoch[which(year < 1840 & year >= 1340)] <- "Maori"
  culturalEpoch[which(year < 1340)] <- "Pre-human"
  return(culturalEpoch)
}

chlgood <- filter(ts,paleoData_variableName == "RABD660670" & paleoData_useInL380National == "TRUE")
allData <- data.frame()
for(i in 1:nrow(chlgood)){
  year <- convertBP2AD(chlgood$age[[i]])
  toplot <- tibble(year = year,
                   depth = chlgood$depth[[i]],
                   RABD660670 = chlgood$paleoData_values[[i]]) %>% 
    mutate(epoch = culturalEpoch(year)) %>% 
    filter(depth > 1, year > 800)
  
  toplot$lake <- chlgood$geo_lakeName[i]
  allData <- bind_rows(allData,toplot)
  
  
  fo <- ggplot(toplot) + 
    geom_line(aes(x = year, y = RABD660670,color = epoch)) +
    scale_color_brewer("Cultural Epoch",palette = "Dark2") + 
    labs(x = "Year AD",y = "RABD 660-670",title = chlgood$geo_siteName[[i]]) +
    theme_bw()+
    theme(legend.position = c(.2,.8))
  
  ggsave(plot = fo,filename = glue::glue("plots/RABD660670/{chlgood$geo_siteName[[i]]}-RABD6606670-CulturalEpoch.pdf"),width = 6,height = 4)
  
  fo2 <- fo + xlim(c(1,2020))
  ggsave(plot = fo2,filename = glue::glue("plots/RABD660670-2k/{chlgood$geo_siteName[[i]]}-RABD6606670-CulturalEpoch.pdf"),width = 6,height = 4)
}


system(glue::glue("rm {file.path(getwd(),'plots','RABD660670','allplots.pdf')}"))
system(glue::glue("cd {file.path(getwd(),'plots','RABD660670')};pdftk *.pdf cat output allplots.pdf"))

  system(glue::glue("rm {file.path(getwd(),'plots','RABD660670-2k','allplots.pdf')}"))
  system(glue::glue("cd {file.path(getwd(),'plots','RABD660670-2k')};pdftk *.pdf cat output allplots.pdf"))

  

  
