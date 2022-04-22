within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Solvent;
function calc_loggamma
  "Calculates decadic log of activity coefficient of water based on Pitzer model"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real loggamma;
protected
  Real log_a = Solvent.calc_log_a(T,p,X);
algorithm
  loggamma :=log_a + log10(IF97.MH2O);
annotation(smoothOrder=20);
end calc_loggamma;
