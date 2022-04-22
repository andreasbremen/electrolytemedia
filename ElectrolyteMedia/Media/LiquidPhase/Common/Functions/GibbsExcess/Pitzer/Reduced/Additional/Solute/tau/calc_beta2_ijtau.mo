within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solute.tau;
function calc_beta2_ijtau "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  output Real[nLifun,nLifun] beta2_ij;
protected
  Real tau = calc_tau(T);
algorithm
  beta2_ij := -Tref/tau^2*(-interactionfun.beta2_ij[:,:,2]*(1/T^2) + interactionfun.beta2_ij[:,:,3]*1/T + interactionfun.beta2_ij[:,:,4] + interactionfun.beta2_ij[:,:,5]*(2*T) + interactionfun.beta2_ij[:,:,6]*(-2/T^3));

end calc_beta2_ijtau;
