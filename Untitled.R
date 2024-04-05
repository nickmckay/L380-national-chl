library(lipdR)

#load in the lipd object
L <- readLipd("~/Download/RichterProtoype.lpd") 

#convert to timeseries object
ts <- as.lipdTsTibble(L)

unique(ts$paleoData_sensorGenus)

this <- dplyr::filter(ts,geo_elevation > 100)




