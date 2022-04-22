within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solvent.pi;
function calc_log_api "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real log_api;
protected
  Real ln_api = calc_ln_api(T,p,X);
algorithm
  log_api :=ln_api*log10(exp(1));

end calc_log_api;
