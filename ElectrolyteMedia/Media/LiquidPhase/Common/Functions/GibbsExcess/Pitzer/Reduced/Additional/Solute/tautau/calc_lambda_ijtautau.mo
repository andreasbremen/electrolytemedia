within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solute.tautau;
function calc_lambda_ijtautau "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  output Real[nLifun,nLifun] lambda_ij;
protected
  Real tau = calc_tau(T);
algorithm
  lambda_ij := 2*Tref/tau^3*(-interactionfun.lambda_ij[:,:,2]*(1/T^2) + interactionfun.lambda_ij[:,:,3]*1/T + interactionfun.lambda_ij[:,:,4] + interactionfun.lambda_ij[:,:,5]*(2*T) + interactionfun.lambda_ij[:,:,6]*(-2/T^3) + interactionfun.lambda_ij[:,:,7]*3*T^2 + interactionfun.lambda_ij[:,:,8]*(-1)/(T-263)^(2))
              + (-Tref/tau^2)^2*(2*interactionfun.lambda_ij[:,:,2]*(1/T^3) - interactionfun.lambda_ij[:,:,3]*1/T^2 + interactionfun.lambda_ij[:,:,5]*2 + interactionfun.lambda_ij[:,:,6]*(6/T^4) + interactionfun.lambda_ij[:,:,7]*6*T + 2*interactionfun.lambda_ij[:,:,8]*1/(T-263)^(3));

end calc_lambda_ijtautau;
