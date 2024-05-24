library(tidyverse)
ts$minAge <- map_dbl(ts$age,min,na.rm = TRUE)
ts$maxYear <- map_dbl(ts$year,max,na.rm = TRUE)
ts$minYear <- map_dbl(ts$year,min,na.rm = TRUE)

unique(ts$inCompilationBeta1_compilationName)


nz <- filter(ts,
             between(geo_latitude,-48, -34),
             between(geo_longitude,166,179),
             minAge < 500 | maxYear > 1450,
             interpretation1_variable == "temperature" |
             interpretation2_variable == "temperature")

nz$year <- map(nz$age,geoChronR::convertBP2AD)

#P <- readLipd("~/Dropbox/lipdverse/Pages2kTemperature/")
pts <- as.lipdTsTibble(P)

nzp <- filter(pts,
             between(geo_latitude,-48, -34),
             between(geo_longitude,166,179),
             interpretation1_variable == "T") |>
  select(-starts_with("has"),-paleoData_inCompilation,-starts_with("pub"))



nzt <- bind_rows(nz,nzp)


nzl <- purrr::transpose(nzt) |>
  tidyTs(age.var = "year") |>
  filter(year > 0)


geoChronR::plotTimeseriesStack(nzl,time.var = "year",color.var = "archiveType")


viol <- nzl |>
  filter(dataSetName %in% c("Aus-Oroko.Cook.2002","Hollywoodcave.Whittaker.2011")) |>
  mutate(era = case_when(
    year < 1250 ~ "PreHuman",
    between(year,1250,1830) ~ "Maori",
    #between(year,1830,1950) ~ "European",
    year > 1830 ~ "European"),
    #era = factor(era,levels = c("PreHuman","Maori","European","Modern"))
    era = factor(era,levels = c("PreHuman","Maori","European"))
  )


ggplot(viol) +
  geom_violin(aes(y = paleoData_values,x = era,fill = dataSetName)) +
  facet_grid(dataSetName ~ .,scales = "free_y") +
  ylab("Proxy Data")



)
