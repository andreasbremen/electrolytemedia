within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Debye.Solvent;
function calc_log_a
  "Calculates decadic logarithm of water activity based on Debye limiting law"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real log_a;
protected
  Real ln_a = calc_ln_a(T,p,X);
algorithm
  log_a :=ln_a*log10(exp(1));

end calc_log_a;
