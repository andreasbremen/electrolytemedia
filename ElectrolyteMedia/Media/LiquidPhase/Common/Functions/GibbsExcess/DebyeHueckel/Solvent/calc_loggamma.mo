within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel.Solvent;
function calc_loggamma
  "calculates decadic log of activity coefficient of water based on extended Debye-Hückel model"
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
