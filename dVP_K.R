# function to calculate the derivative of saturation vapor
# pressure at Ta (Flatau et al., 1992)
require(pracma)
dVP_K=function(Ta){
    b = c(0.381294516e-16, -0.114596802e-13, -0.788037859e-12,
	 0.404125005e-09, 0.103354611e-06, 0.121211669e-04,
	 0.794683137e-03, 0.286064092e-01, 0.444017302)
    To=273.16
    x = Ta-To 
    z = polyval(b,x)
    z =z*100 #saturating vapor pressue air [Pa] /100 gives mbar.
    z      
}