within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Solute;
function calc_psi_ijk
  "Calculates temperature dependent Pitzer polynomial for triple charged interaction"
  input Modelica.SIunits.Temperature T;
  output Real[nLifun,nLifun,nLifun] psi_ijk;
algorithm
  psi_ijk := interactionfun.psi_ijk[:,:,:,1] + interactionfun.psi_ijk[:,:,:,2]*(1/T-1/Tref) + interactionfun.psi_ijk[:,:,:,3]*log(T/Tref) + interactionfun.psi_ijk[:,:,:,4]*(T-Tref) + interactionfun.psi_ijk[:,:,:,5]*(T^2-Tref^2) + interactionfun.psi_ijk[:,:,:,6]*(1/T^2-1/Tref^2);

end calc_psi_ijk;
