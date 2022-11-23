
calcFlux <- function(paleoData_values,bulkDensity, sedRate){
   calibration <-  calcConc(paleoData_values)
 
  return(calibration * bulkDensity * sedRate)
}

calcConc <- function(paleoData_values){
  ten_root <- function(x) x ^ (1/10)
  ten <- function(x) x ^ 10
  minus_root <- function(x) (x-.98) ^ (1/10)
  minus <- function(x) (x^10)+.98
  minRoot <- minus_root(paleoData_values)
  
  #calibration <-  ten(-0.4707 +(2.3629 * minRoot)) Ethan's
  calibration <-  ten(-0.8238 +(2.7233 * minRoot)) #Nick retry
  
  return(calibration)
}

estBulkDensity <- function(dblf){
  return(log(dblf+1)/12+.9)
}

culturalEpoch <- function(year){
  culturalEpoch <- matrix(NA,nrow = length(year))
  culturalEpoch[which(year >= 1840)] <- "European"
  culturalEpoch[which(year < 1840 & year >= 1340)] <- "Maori"
  culturalEpoch[which(year < 1340)] <- "Pre-human"
  culturalEpoch <- factor(culturalEpoch,levels = c("Pre-human","Maori","European"))
  return(culturalEpoch)
}

deltas <- function(vec){
  return(c(diff(vec),NA))
}

whichIsland <- function(lat,lon){
 if(lon > 174 | lat > -40){
   return("North Island")
 }else{
   return("South Island")
 }
}



