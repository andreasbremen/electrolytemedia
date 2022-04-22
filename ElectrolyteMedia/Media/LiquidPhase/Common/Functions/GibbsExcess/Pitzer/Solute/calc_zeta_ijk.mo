within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Solute;
function calc_zeta_ijk
  "Calculates temperature dependent Pitzer polynomial for neutral-charged-charged interaction"
  input Modelica.SIunits.Temperature T;
  output Real[nLifun,nLifun,nLifun] zeta_ijk;
algorithm
  zeta_ijk := interactionfun.zeta_ijk[:,:,:,1] + interactionfun.zeta_ijk[:,:,:,2]*(1/T-1/Tref) + interactionfun.zeta_ijk[:,:,:,3]*log(T/Tref) + interactionfun.zeta_ijk[:,:,:,4]*(T-Tref) + interactionfun.zeta_ijk[:,:,:,5]*(T^2-Tref^2) + interactionfun.zeta_ijk[:,:,:,6]*(1/T^2-1/Tref^2) + interactionfun.zeta_ijk[:,:,:,7]*T^3 + interactionfun.zeta_ijk[:,:,:,8]*1/(T-263);

end calc_zeta_ijk;
