while (abs(dc)>=1e-07){ #loop while the absolute value of dc > or = 0.1 ppm
    
    niter=niter+1 # add to counter
    
    cc0=cc0-dc/10 # update the value of cc0 based on dc from the previous loop
    
    pCO2c=cc0*ptot # convert cc0 to partial pressure used in script
    
    source('doPS_pCO2c.R')  # run script
    
    gbc = gb /1.4
    gsc = g0 /1.6
    gtc = 1/(1/gsc+1/gbc)
    cc1=ca - assimn/gtc # estimate mole fraction CO2 inside the leaf from the diffusion equation

  dc=cc0-cc1 #calculates the error 
		# (difference between the guess cc0 and the calculated value)
  #make new estimate of cc0
  ga=assimn/cc0
  cc0=ca * gtc /(gtc+ga)
    if(niter>200){
      break
    } # added to avoid wavering between two possible values
}

