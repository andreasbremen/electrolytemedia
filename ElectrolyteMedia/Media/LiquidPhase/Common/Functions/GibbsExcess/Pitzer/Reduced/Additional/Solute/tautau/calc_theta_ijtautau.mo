within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solute.tautau;
function calc_theta_ijtautau "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  output Real[nLifun,nLifun] theta_ij;
protected
  Real tau = calc_tau(T);
algorithm
  theta_ij := 2*Tref/tau^3*(-interactionfun.theta_ij[:,:,2]*(1/T^2) + interactionfun.theta_ij[:,:,3]*1/T + interactionfun.theta_ij[:,:,4] + interactionfun.theta_ij[:,:,5]*(2*T) + interactionfun.theta_ij[:,:,6]*(-2/T^3))
              + (-Tref/tau^2)^2*(2*interactionfun.theta_ij[:,:,2]*(1/T^3) - interactionfun.theta_ij[:,:,3]*1/T^2 + interactionfun.theta_ij[:,:,5]*2 + interactionfun.theta_ij[:,:,6]*(6/T^4));

end calc_theta_ijtautau;
