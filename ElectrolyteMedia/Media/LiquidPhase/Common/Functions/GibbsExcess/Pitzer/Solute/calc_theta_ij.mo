within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Solute;
function calc_theta_ij
  "Calculates temperature dependent Pitzer polynomial for binary anion-anion/cation-cation interaction"
  input Modelica.SIunits.Temperature T;
  output Real[nLifun,nLifun] theta_ij;
algorithm
  theta_ij := interactionfun.theta_ij[:,:,1] + interactionfun.theta_ij[:,:,2]*(1/T-1/Tref) + interactionfun.theta_ij[:,:,3]*log(T/Tref) + interactionfun.theta_ij[:,:,4]*(T-Tref) + interactionfun.theta_ij[:,:,5]*(T^2-Tref^2) + interactionfun.theta_ij[:,:,6]*(1/T^2-1/Tref^2);

end calc_theta_ij;
