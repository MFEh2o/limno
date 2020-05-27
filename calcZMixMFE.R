#Function to calculate thermocline depth from manual temp profiles in MFE database
#CTS 3 Oct 2017

#Based partly on calcZMixDens.r, a function I wrote several years ago for analysis of
#data for Solomon et al. 2013 L&O

calcZMixMFE <- function(d,thresh=0.1) {
  
  #This function converts temperature profile to density and determines the thermocline
  #as the shallowest depth at which the rate of change of density with depth exceeds a threshold.
  
  #Arguments:
  # d         Input data.frame with (minimally) the following columns:
  #             -temp - water temperature in C
  #             -depthTop - depth (m) at which temperature was measured
  # thresh    Threshold density change to indicate end of mixed layer. Units (kg m-3) m-1.
  #           Defaults to 0.1 (default in rLakeAnalyzer thermo.depth function)
  
  #Value:
  # zMix      The depth (m) of the thermocline
  
  
  #Convert temp data to density
  #Density of water (kg m-3) as function of temp from McCutcheon (1999)
  #Note there is a different formula if salinity is appreciable; formula below ignores that
  d$dens <- 1000*(1-((d$temp+288.9414)/(508929.2*(d$temp+68.12963)))*(d$temp-3.9863)^2)
  
  #If there are not at least two non-NA dens values, return NA for zMix
  if (length(which(is.na(d$dens)==FALSE))<2) {
    zMix <- NA
  } else {
    
    #Linearly interpolate the density profile (50 points between min(depthTop) and max(depthTop)
    interpProfile <- approx(x=d$depthTop,y=d$dens)
    
    #Calculate change in density with depth
    dDdZ <- diff(interpProfile$y)/diff(interpProfile$x)
    
    #If no rate of change exceeds thresh, return max depth for zMix
    if (all(dDdZ<thresh)) {
      zMix <- max(d$depthTop)
    } else {
      
      #Otherwise, determine the shallowest depth at which dDdZ exceeds threshold density difference
      whichMin <- min(which(dDdZ>thresh))
      zMix <- mean(c(interpProfile$x[whichMin],interpProfile$x[(whichMin+1)]))
    }
      
  }
  
  #Return zMix
  zMix <- as.data.frame(zMix)
  return(zMix)
}
  