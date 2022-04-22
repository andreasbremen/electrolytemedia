within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Solute;
function calc_beta2_ij
  "Calculates temperature dependent Pitzer polynomial for binary cation-anion interaction"
  input Modelica.SIunits.Temperature T;
  output Real[nLifun,nLifun] beta2_ij;
algorithm
  beta2_ij := interactionfun.beta2_ij[:,:,1] + interactionfun.beta2_ij[:,:,2]*(1/T-1/Tref) + interactionfun.beta2_ij[:,:,3]*log(T/Tref) + interactionfun.beta2_ij[:,:,4]*(T-Tref) + interactionfun.beta2_ij[:,:,5]*(T^2-Tref^2) + interactionfun.beta2_ij[:,:,6]*(1/T^2-1/Tref^2);

end calc_beta2_ij;
