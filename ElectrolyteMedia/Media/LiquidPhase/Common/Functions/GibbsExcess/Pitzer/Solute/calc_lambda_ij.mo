within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Solute;
function calc_lambda_ij
  "Calculates temperature dependent Pitzer polynomial for neutral-(charged ion/neutral) interaction"
  input Modelica.SIunits.Temperature T;
  output Real[nLifun,nLifun] lambda_ij;
algorithm
  lambda_ij := interactionfun.lambda_ij[:,:,1] + interactionfun.lambda_ij[:,:,2]*(1/T-1/Tref) + interactionfun.lambda_ij[:,:,3]*log(T/Tref) + interactionfun.lambda_ij[:,:,4]*(T-Tref) + interactionfun.lambda_ij[:,:,5]*(T^2-Tref^2) + interactionfun.lambda_ij[:,:,6]*(1/T^2-1/Tref^2) + interactionfun.lambda_ij[:,:,7]*T^3 + interactionfun.lambda_ij[:,:,8]*1/(T-263);
end calc_lambda_ij;
