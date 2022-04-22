within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solute.tautau;
function calc_beta2_ijtautau "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  output Real[nLifun,nLifun] beta2_ij;
protected
  Real tau = calc_tau(T);
algorithm
  beta2_ij := 2*Tref/tau^3*(-interactionfun.beta2_ij[:,:,2]*(1/T^2) + interactionfun.beta2_ij[:,:,3]*1/T + interactionfun.beta2_ij[:,:,4] + interactionfun.beta2_ij[:,:,5]*(2*T) + interactionfun.beta2_ij[:,:,6]*(-2/T^3))
              + (-Tref/tau^2)^2*(2*interactionfun.beta2_ij[:,:,2]*(1/T^3) - interactionfun.beta2_ij[:,:,3]*1/T^2 + interactionfun.beta2_ij[:,:,5]*2 + interactionfun.beta2_ij[:,:,6]*(6/T^4));

end calc_beta2_ijtautau;
