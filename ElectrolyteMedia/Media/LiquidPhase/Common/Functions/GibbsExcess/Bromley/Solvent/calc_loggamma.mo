within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Bromley.Solvent;
function calc_loggamma
  "calculates decadic log of activity coefficient of water based on Bromleys method with extension to multielectrolyte systems by Meissner and Kusik"

  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nLfun] X;

  output Real loggamma;
protected
  Real log_a=calc_log_a(T,p,X);
algorithm
  loggamma :=log_a + log10(IF97.MH2O);
annotation(smoothOrder=20);
end calc_loggamma;
