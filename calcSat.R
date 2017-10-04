#Function to calculate saturation concentration of oxygen dissolved in water
#CTS 4 Oct 2017 based on code in the metabFunc_v6.r script used for Solomon et al. 2013 L&O analysis

calcSat <- function(temp,atmPres) {
  
  #This function calculates DO saturation following Weiss 1970 Deep Sea Res. 17:721-735; 
  #simplified since salinity=0.
  #Equation is ln DO = A1 + A2 100/T + A3 ln T/100 + A4 T/100
  
  #Arguments:
  # temp    Water temperature (C)
  # atmPres Atmospheric pressure (mm Hg); can be measured or a local average given elevation
  
  #Value:
  # DOSat   Cocentration of dissolved oxygen (mg L-1) at saturation given temperature and pressure
  
  #Convert water temp to Kelvin
  tempK <- temp + 273.15
  
  #Weiss equation
  A1 <- -173.4292;  A2 <- 249.6339;  A3 <- 143.3483;  A4 <- -21.8492
  DOSat <- exp(((A1 + (A2*100/tempK) + A3*log(tempK/100) + A4*(tempK/100))))
  
  #Correction for local average atmospheric pressure
  u <- 10^(8.10765 - (1750.286/(235+temp)))
  DOSat <- (DOSat*((atmPres-u)/(760-u)))   #ml/L
  DOSat <- DOSat/1000                      #L/L
  
  #Convert using standard temperature and pressure. 
  #Similar to calculating saturation DO at STP in ml/L, converting to mg?L (at STP),
  #and then doing the above temperature and pressure conversions.
  R <- 0.082057  #L atm deg-1 mol-1
  O2molWt <- 15.999*2
  convFactor <- O2molWt*(1/R)*(1/273.15)*(760/760) #g/L
  DOSat <- DOSat*convFactor*1000                   #mg/L
  
  #Return DOSat
  return(DOSat)
}