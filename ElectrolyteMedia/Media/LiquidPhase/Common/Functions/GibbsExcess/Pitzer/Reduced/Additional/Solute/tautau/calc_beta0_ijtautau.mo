within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solute.tautau;
function calc_beta0_ijtautau "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  output Real[nLifun,nLifun] beta0_ij;
protected
  Real tau = calc_tau(T);
algorithm
  beta0_ij := 2*Tref/tau^3*(-interactionfun.beta0_ij[:,:,2]*(1/T^2) + interactionfun.beta0_ij[:,:,3]*1/T + interactionfun.beta0_ij[:,:,4] + interactionfun.beta0_ij[:,:,5]*(2*T) + interactionfun.beta0_ij[:,:,6]*(-2/T^3) + interactionfun.beta0_ij[:,:,7]*3*T^2 + interactionfun.beta0_ij[:,:,8]*(-1)/(T-263)^(2))
              + (-Tref/tau^2)^2*(2*interactionfun.beta0_ij[:,:,2]*(1/T^3) - interactionfun.beta0_ij[:,:,3]*1/T^2 + interactionfun.beta0_ij[:,:,5]*2 + interactionfun.beta0_ij[:,:,6]*(6/T^4) + interactionfun.beta0_ij[:,:,7]*6*T + 2*interactionfun.beta0_ij[:,:,8]*1/(T-263)^(3));

end calc_beta0_ijtautau;
