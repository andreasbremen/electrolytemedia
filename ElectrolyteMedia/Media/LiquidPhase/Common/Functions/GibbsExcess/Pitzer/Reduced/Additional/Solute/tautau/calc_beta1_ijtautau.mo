within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solute.tautau;
function calc_beta1_ijtautau "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  output Real[nLifun,nLifun] beta1_ij;
protected
  Real tau = calc_tau(T);
algorithm
  beta1_ij := 2*Tref/tau^3*(-interactionfun.beta1_ij[:,:,2]*(1/T^2) + interactionfun.beta1_ij[:,:,3]*1/T + interactionfun.beta1_ij[:,:,4] + interactionfun.beta1_ij[:,:,5]*(2*T) + interactionfun.beta1_ij[:,:,6]*(-2/T^3) + interactionfun.beta1_ij[:,:,7]*3*T^2 + interactionfun.beta1_ij[:,:,8]*(-1)/(T-263)^(2))
              + (-Tref/tau^2)^2*(2*interactionfun.beta1_ij[:,:,2]*(1/T^3) - interactionfun.beta1_ij[:,:,3]*1/T^2 + interactionfun.beta1_ij[:,:,5]*2 + interactionfun.beta1_ij[:,:,6]*(6/T^4) + interactionfun.beta1_ij[:,:,7]*6*T + 2*interactionfun.beta1_ij[:,:,8]*1/(T-263)^(3));

end calc_beta1_ijtautau;
