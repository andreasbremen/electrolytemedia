within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel.Solvent;
function calc_log_a
  "calculates decadic logarithm of water activity based on extended Debye Hückel model"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real log_a;
protected
  Real ln_a = calc_ln_a(T,p,X);
algorithm
  log_a :=ln_a*log10(exp(1));

end calc_log_a;
