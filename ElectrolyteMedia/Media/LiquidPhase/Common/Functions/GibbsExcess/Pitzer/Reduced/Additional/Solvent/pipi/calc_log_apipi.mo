within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solvent.pipi;
function calc_log_apipi "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real log_apipi;
protected
  Real ln_apipi = calc_ln_apipi(T,p,X);
algorithm
  log_apipi :=ln_apipi*log10(exp(1));

end calc_log_apipi;
