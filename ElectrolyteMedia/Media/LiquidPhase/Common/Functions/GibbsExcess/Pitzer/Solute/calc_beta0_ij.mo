within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Solute;
function calc_beta0_ij
  "Calculates temperature dependent Pitzer polynomial for binary cation-anion interaction"
  input Modelica.SIunits.Temperature T;
  output Real[nLifun,nLifun] beta0_ij;
algorithm
  beta0_ij := interactionfun.beta0_ij[:,:,1] + interactionfun.beta0_ij[:,:,2]*(1/T-1/Tref) + interactionfun.beta0_ij[:,:,3]*log(T/Tref) + interactionfun.beta0_ij[:,:,4]*(T-Tref) + interactionfun.beta0_ij[:,:,5]*(T^2-Tref^2) + interactionfun.beta0_ij[:,:,6]*(1/T^2-1/Tref^2) + interactionfun.beta0_ij[:,:,7]*T^3 + interactionfun.beta0_ij[:,:,8]*1/(T-263);

end calc_beta0_ij;
