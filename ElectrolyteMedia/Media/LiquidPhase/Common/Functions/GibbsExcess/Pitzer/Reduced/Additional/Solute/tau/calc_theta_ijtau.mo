within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solute.tau;
function calc_theta_ijtau "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  output Real[nLifun,nLifun] theta_ij;
protected
  Real tau = calc_tau(T);
algorithm
  theta_ij := -Tref/tau^2*(-interactionfun.theta_ij[:,:,2]*(1/T^2) + interactionfun.theta_ij[:,:,3]*1/T + interactionfun.theta_ij[:,:,4] + interactionfun.theta_ij[:,:,5]*(2*T) + interactionfun.theta_ij[:,:,6]*(-2/T^3));

end calc_theta_ijtau;
