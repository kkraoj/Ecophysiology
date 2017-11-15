#function to calculate the saturation vapor pressure of water at Ta 
#	(Flatau et al., 1992)
require(pracma) # evaluate a polynomial with coefficients of order n to 0
VP_K=function(Ta){
    a = c(0.209339997e-15, -0.373208410e-12, 0.892344772e-10,
	0.196237241e-07, 0.305903558e-05, 0.264461437e-03, 0.143064234e-01,
	0.444007856, 6.11213476)
    To=273.16
    x = (Ta-To) 
    z = polyval(a,x)   #saturating vapor pressue air [Pa]
    z =z*100 		#divide by 100 to get mbar
    z
}     
                        
