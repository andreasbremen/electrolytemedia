within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Debye.Solvent;
function calc_loggamma
  "Calculates decadic log of activity coefficient of water based on Debye limiting law"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real loggamma;
protected
  Real log_a = calc_log_a(T,p,X);
algorithm
  loggamma :=log_a + log10(IF97.MH2O);
annotation(smoothOrder=20);
end calc_loggamma;
