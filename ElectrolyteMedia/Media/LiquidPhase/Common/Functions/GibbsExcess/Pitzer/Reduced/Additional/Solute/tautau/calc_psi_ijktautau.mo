within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solute.tautau;
function calc_psi_ijktautau "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  output Real[nLifun,nLifun,nLifun] psi_ijk;
protected
  Real tau = calc_tau(T);
algorithm
  psi_ijk := 2*Tref/tau^3*(-interactionfun.psi_ijk[:,:,:,2]*(1/T^2) + interactionfun.psi_ijk[:,:,:,3]*1/T + interactionfun.psi_ijk[:,:,:,4] + interactionfun.psi_ijk[:,:,:,5]*(2*T) + interactionfun.psi_ijk[:,:,:,6]*(-2/T^3))
              + (-Tref/tau^2)^2*(2*interactionfun.psi_ijk[:,:,:,2]*(1/T^3) - interactionfun.psi_ijk[:,:,:,3]*1/T^2 + interactionfun.psi_ijk[:,:,:,5]*2 + interactionfun.psi_ijk[:,:,:,6]*(6/T^4));

end calc_psi_ijktautau;
