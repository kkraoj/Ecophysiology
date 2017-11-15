
while (abs(Dt)>=0.01){
  nitEB=nitEB+1
    #  t1 = new guess for leaf temperature
    t1=t1+Dt
    #  Re = long wave radiation emitted (W m^-2) top + bottom of leaf
    tk = t1 + 273.16
    Re = 2 * ems * sbc * tk^4
    #  Rn = net radiation (Ra - Re)
    Rn = Ra - Re
    # derivative of radiative dissipation two sides of leaf
    dR = 8 * ems * sbc * tk^3
    #  H = sensible heat transfer to air(W m^-2)
    gh=0.91 * gb
    H = cpa * (t1 - ta) * 2 * gh # two sides
    #  dH = derivative of H
    dH = cpa * 2 * gh
    #  E = evaporation rate (mol m^-2 s^-1)
    esat = VP_K(tk)
    E = gv *(esat - ea)/ptot
    #  LE = latent heat (W m^-2)
    LE = E * lhv
    #  dle = derivative of lE
    desat=dVP_K(tk)
    dLE = lhv * gv * desat/ptot
    #  Y = error in energy balance calculation
    Y = Rn - H - LE
    #  t1 = guess for leaf temperature
    #  Dt = Delta-temperature
    Dt = Y/(dR + dH + dLE)#  left off here to show progress of iteration
  if (nitEB > 100){
    stopped=1
    break
  } # workaround to keep the loop from spiraling out of control
        #   as it will do otherwise for certain parameter combinations
}
