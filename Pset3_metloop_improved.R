#### DoFlux script, modified to loop over the Manaus met dataset

met = read.csv('Manausmet.csv')
head(met)

# Make the date into a time object: not necessary for this assignment, but useful
# R has several ways to deal with time objects 
library(lubridate)
met$date=as.POSIXct(as.character(met$nymd),format='%y%m%d%H',tz='America/Manaus') 
met$hour=hour(met$date)
met$day=day(met$date)


# Required package to do Matlab-style evaluations
library(pracma)

# Get functions you will need later
source('VP_K.R')
source('dVP_K.R')

# Model configuration inputs
rstfac=1
# cif=1.5
# ty=1 
ca=365e-06
tks=299.8

# Global parameters
Ps_params=read.table('Ps_params_header.csv',sep=',',header=T)


# ----------------------------------------------------------
#   Constants

#  alb = reflectance of leaf to SW radiation 
alb = 0.17

#  ems = longwave emissivity
ems = 0.95

#  sbc = Stefan-Boltzman constant
sbc = 5.6703e-8

#  lhv = latent heat of vaporization (joules mol^-1)
lhv = 44e3

#  cpa = molar specific heat of air (joules mol^-1 oK^-1)   
cpa = 29.3

#  pc = psychrometrer constant (oK^-1) (used in Penman-Monteith)
pc=29.3/lhv


  c4=Ps_params[ty,1]   		#Photosynthetic type c3=0 c4=1
  vmax_umol=Ps_params[ty,2] * cif   #Rubisco, NOTE: mol m-2 s-1 
  smax=Ps_params[ty,3]        #C3 POS factor (unused in C4)
  zkc_p=Ps_params[ty,4]       #C3 Rubisco CO2 Km (carboxylation)      
  zko_p=Ps_params[ty,5]       #C3 Rubisco O2 Km (oxygenation) 
  rrkk_p=Ps_params[ty,6]      #C4 rrkk parameter (unused in C3)    
  leaf_Tc =Ps_params[ty,7]    # default Input leaf temperature,K
  #apar_umol=Ps_params[ty,8]   # default input par (umo m-2 s-2)   
  #Ca_ppm=Ps_params[ty,9]      # default atmospheric CO2 (ppm) 
  ptot=Ps_params[ty,10]       # total atmospheric pressure    
  pO2m=Ps_params[ty,11]       # partial pressure of O2    
  PARalb=Ps_params[ty,12]     # leaf or canopy albedo in PAR band 
  emis=Ps_params[ty,13]       # longwave emissivity   
  nitratefac=Ps_params[ty,14] # nitrogen nutrition scaling factor unused  
  rstfac1=Ps_params[ty,15]    # Water stress factor (0-1) 
  effcon=Ps_params[ty,16]     # Intrinsic quantum efficiency  
  btheta=Ps_params[ty,17]     # Curvature parameter used in quadratic solution    
  hhti=Ps_params[ty,18]       # High temperature inhibition factor, PS    
  hlti=Ps_params[ty,19]       # Low temperature inhibition factor, PS     
  respcp=Ps_params[ty,20]     # Factor relating respiration rate to Rubisco vmax  
  atheta=Ps_params[ty,21]     # Curvature Parameter used in second quadratic  
  shti=Ps_params[ty,22]       # Parameter used in temperature response of PS (normally unchanged)   
  slti=Ps_params[ty,23]       # Parameter used in temperature response, PS (normally unchanged)   
  trda=Ps_params[ty,24]       # Temperature response of respiration   
  trdm=Ps_params[ty,25]       # Temperature response of respiration   
  gradm=Ps_params[ty,26]      # Stomata slope parameter   
  bintc=Ps_params[ty,27]      # Stomata intercept parameter
  
met$A=rep(NA)
met$Rn=rep(NA)
met$H=rep(NA)
met$LE=rep(NA)

for (i in seq(length(met$nymd))){
  swdown=met$swdown[i]
  zlwd=met$zlwd[i]
  em=met$em[i]
  tm=met$tm[i]
  um=met$um[i]

  #convert and/or calculate needed time step inputs
  ea= em * 100 
  ta=tm -273.15
  ts=tks-273.15
  SWd=swdown
  LWd=zlwd
  gb=0.147* sqrt(um/0.1)
  
  Ra = SWd * (1 - alb) + LWd * ems + ems * sbc *(ts+273.16)^4
  # ----------------------------------------------------------
  # calculate apar
  apar=(1-PARalb) * SWd *0.5 /2.35e05 
  #-----------------------------------------------------------
  # start solution by guessing initial looping variables and escape criteria
  #-----------------------------------------------------------
  t0 = ta
  Dt=0.1
  cc0 = (0.8 - 0.4 * c4) * ca
  dc=0.1e-06
  g0=0.5 # mol m-2s-1)
  dg=0.01
  niter=0
  nitEB=0 # separate counter for energy budget
  while (abs(dg)>=.001){
    niter=niter+1
    tleaf=t0
    dc=1e-06
    #calculate assimn as f(g0)
    source('doPS_g0.R')
    #estimate a new conductance based on BB eqn.
    esat=VP_K(t0+273.16)
    rv=1/g0+1/gb
    gv=1/rv
    E=gv * (esat-ea)/ptot
    es=ea+ptot * E /gb
    hs=es/esat
    cs = ca - 1.4 * assimn / gb
    g1 = gradm * assimn * hs / cs+bintc
    dg=g0-g1
    #run the energy budget script to get a new t0 compatible with the new g1
    t1=t0
    Dt=0.01
    source('doEB_g1.R')
    
    if (Y > 5 || assimn < -10^-5){
      assimn=NA
      Rn=NA
      H=NA
      LE=NA
      break
    } 
    # Stop the loop and return NA if values make no sense
    t0=t1
    g0=g1
    # Send back to begining with new leaf temperature and conductance
  }

  met$A[i]=assimn*10^6 # convert to umol to make the plots nicer
  met$Rn[i]=Rn
  met$H[i]=H
  met$LE[i]=LE
}

