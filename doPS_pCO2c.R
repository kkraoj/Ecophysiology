
# Photosynthesis (assimn) at the leaf cell scale - using pCO2c in Pascals.

# The following unputs are required from the calling script:

#
# 1) The data table "PS_params". Each row contains a complete set of parameters

# for a SiB biome type. Any custom paramaterizations can be added by making

# additional rows using the comma-separated spread sheet PS_params.csv.

#
# 2) the variable "ty"  used to select the desired row.

#
# 3) the variable "pCO2c" in Pa. The is the partial pressure of CO2 inside

#      the chloroplast.

#
# 4) the variable "apar"  in mol m-2 s-1. This is the is the absorbed flux of 

#      photosyntheticaly active radiation.

#
# 5) the variable "tleaf"  in oC. This is the leaf temperature.

#tleaf=leaf_Tc


# 
# This script calculates the rate of net CO2 assimilation "asimin" (mol m-2 s-1)

#------------------------------------------------------------------------
  

  
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  
# Biome dependent physiological properties (SE-96,Table 5) for
  
# biome 1: broadleaf-evergreen trees (a.k.a. tropical broadleaf)
  
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
            
       
#un-comment if you want to enter parameter file here  

#     c4=Ps_params[ty,1]   		
#Photosynthetic type c3=0 c4=1
#     vmax_umol=Ps_params[ty,2]   
#Rubisco, NOTE: mol m-2 s-1 
#     smax=Ps_params[ty,3]        
#C3 POS factor (unused in C4)
#     zkc_p=Ps_params[ty,4]       
#C3 Rubisco CO2 Km (carboxylation)      

#     zko_p=Ps_params[ty,5]       
#C3 Rubisco O2 Km (oxygenation) 

#     rrkk_p=Ps_params[ty,6]      
#C4 rrkk parameter (unused in C3)    

#     leaf_Tc =Ps_params[ty,7]    
# default Input leaf temperature,oK
#     apar_umol=Ps_params[ty,8]   
# default input par (umo m-2 s-2)   

#     #Ca_ppm=Ps_params[ty,9]      
# default atmospheric CO2 (ppm) 
#     ptot=Ps_params[ty,10]       
# total atmospheric pressure    
#     pO2m=Ps_params[ty,11]       
# partial pressure of O2    
#     PARalb=Ps_params[ty,12]     
# leaf or canopy albedo in PAR band 
#     emis=Ps_params[ty,13]       
# longwave emissivity   
#     nitratefac=Ps_params[ty,14] 
# nitrogen nutrition scaling factor unused  

#     rstfac1=Ps_params[ty,15]    
# Water stress factor (0-1) 

#     effcon=Ps_params[ty,16]     
# Intrinsic quantum efficiency  

#     btheta=Ps_params[ty,17]     
# Curvature parameter used in quadratic solution    

#     hhti=Ps_params[ty,18]       
# High temperature inhibition factor, PS    

#     hlti=Ps_params[ty,19]       
# Low temperature inhibition factor, PS     

#     respcp=Ps_params[ty,20]     
# Factor relating respiration rate to Rubisco vmax  

#     atheta=Ps_params[ty,21]     
# Curvature Parameter used in second quadratic  

#     shti=Ps_params[ty,22]       
# Parameter used in temperature response of PS (normally unchanged)   
#     slti=Ps_params[ty,23]       # Parameter used in temperature response, PS (normally unchanged)   
#     trda=Ps_params[ty,24]       # Temperature response of respiration   
#     trdm=Ps_params[ty,25]       # Temperature response of respiration   
#     gradm=Ps_params[ty,26]      # Stomata slope parameter   
#     bintc=Ps_params[ty,27]      # Stomata intercept parameter
    
 
  #------------------------------------------------------------------
  # Make any temporary changes in the above parameter values by editing
  # here.  For example.
  # 
  # vmax_umol = 80 
  #
  #  would change the value of this parameter from that read in above
  # to this new value.
  #------------------------------------------------------------------
  
  #------------------------------------------------------------------
  # Don't change anything below here
  #------------------------------------------------------------------
  
  #------------------------------------------------------------------
  
# Setup common to both C3 and C4 calculations
  #------------------------------------------------------------------
  
  # Calculate temperature and stress effects
  # See SE-96A, pg 701, (C17)
  # The reference temperature is 25oC for all kinetic constants
  # See also SvC, Biochemical Models of Leaf Photosynthesis, Table 2.3
  
  tk=tleaf+273.16 #leaf temperature to kelvin
  
  # Stress limits on Vmax and Respiration
  templ = 1 + exp(slti * (hlti-tk)) # low temperature, C4 only
  temph = 1 + exp(shti * (tk-hhti)) # high temperature, C3 & C4
  rstfac4 = 1/(templ * temph)       # unused with C3

  # Exponent of the Q10 function (Q10 temperature coefficient)
  #   Sellers et al '96, Part 2, Table 5 (unitless)
  qt = (tk - 298) / 10
  
  #  Use Vmax Rubisco from the input Ps_prarams change to mol m-2 s-1
  vmax=vmax_umol*1e-06
  
  #------------------------------------------------------------------
  # Use c4 parameter to select c3 or c4 calculation
  #------------------------------------------------------------------
  
if (c4==0){  
  # Begin C3 Photosynthesis Calculation
  
  # Specificity factor of Rubisco at tk (tau, Collatz,91 A3)
  #   Sellers et al '96, Part 2, Table 5 (unitless)
  spfy = 2600 * (0.57^qt)

  # Km(carboxylation, CO2) of Rubisco at this temperature
  #   Sellers et al '96, Part 2, Table 5 (30 Pa.)
  #   SvC Biochemical Models book, Table 2.3 (260 ubar = 26 Pa)
  zkc = zkc_p * (2.1^qt)

  # Km(oxygenation = photorespiratiot, O2) of Rubisco at tk
  #   Sellers et al '96, Part 2, Table 5 (30000 Pa.)
  #   SvC Biochemical Models book, Table 2.3 (179 mbar = 17900 Pa)
  zko = zko_p * (1.2^qt)

  # Rubisco Vmax at this temperature  
  vm = vmax * (2.1^qt)
  
  # Adjust vm for temperature & water stress
  # rstfac1 is water stress parameter (see phosib.F for clue)
  vm = vm * rstfac1 / temph

  # Scale respiration to vmax
  respc = respcp * vmax * (1.8^qt)
  
  # High temperature limit on respiration
  respc = respc / (1 + exp(trda * (tk - trdm)))

  # compensation point, "gamma star" (Collatz91 A3)
  gammas = 0.5 * pO2m / spfy
  gstar=gammas /ptot # gammas in mol fraction 
  

  # oxygen effect on zkc
  rrkk = zkc * (1.0 + pO2m / zko)
  
  # initial slope of A/Ci curve at gammas vC#F 
  daci=ptot * vm /(gammas+rrkk)

  # estimate mesophyll (cell wall) conductance 
  #gmeso = vmax * 1e4

  # actual CO2 compensation point (can be used in iteration)
  
   gamma=gammas/ptot+respc / daci
  
  # calculate limiting rates Collatz'90 A5, A2, & A7
  
  omc = vm * (pCO2c - gammas) / (pCO2c + rrkk) 
  ome = apar * effcon * (pCO2c - gammas) / (pCO2c + 2 * gammas)           
  oms = vmax * rstfac * 1.8 ^qt / 2  
}else{
   
  #Begin C4 calculation after Collatz '92
    
  # Rubisco Vmax at this temperature  
  vm = vmax * (2.1^qt)

  # Adjust vm for temperature & water stress
  # rstfac1 is water stress parameter (see phosib.F for clue)
  vm = vm * rstfac1 * rstfac4

  # Scale respiration to vmax
  respc = respcp * vmax * (1.8^qt)

  # High temperature limit on respiration
  respc = respc / (1.0 + exp(trda * (tk - trdm)))

  # first order rate constant for PEP case
  # the rrkk_parameter affects the initial slope of the C02xA curve 
  # 
  rrkk = 1.8^qt * vmax/5
  daci=ptot * rrkk

  # estimate mesophyll (cell wall) conductance (not used here)
  gmeso = vmax * 3e4

  #compensation point
  gammas=0.1 #assmm 1ppm for C4
  gstar=gammas /ptot # gammas in mol fraction 

  # Calculate carbon assimilation limiting rates
  
  omc = vm
  ome = apar * effcon
  oms = rrkk * pCO2c
} 


#------------------------------------------------------------------
# Solve for assim using omc, ome and oms calculated above - common to both
# c3 and c4 calculations
#------------------------------------------------------------------
  # jump out of solution at or below the light compensation point
  
  if (ome <= respc){
      assim=ome
}else{
        #------------------------------------------------------------------  
        #    Solve the first quadratic for omp
        #    atheta*omp^2 - (omc + ome)*omp + ome*omc=0
        #------------------------------------------------------------------
        

        #a = c( omc*ome, -(omc + ome), atheta )
        #a1 = Re(polyroot(a)) # base R method of solving a polynomial
	  # R lists polynomial terms in order of increasing degree (0th to nth)
	  	# MATLAB-style, with the "Practical Numerical Math" library: 
		  require(pracma) 
	  	  a = c( atheta, -(omc + ome), omc*ome ) 
	  	# ordered nth to 0th can also use rev()
	  	  a1 = roots(a)

        if (pCO2c > gammas){
            omp = min(a1)
	  }else{
            omp = max(a1)
        }
        
        #------------------------------------------------------------------
        #       Solve the second quadratic for assim
        #       btheta*assim^2 - (omp+oms)*assim + omp*oms = 0
        #------------------------------------------------------------------                                         
        
        b =c(omp*oms,-(omp + oms),btheta)
        b1 = Re(polyroot(b))
        assim = min(b1)
           
        
  }
  assimn =  (assim - respc)
  assimn
  
#------------------------------------------------------------------------


