within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solvent.pitau;
function calc_loggammapitau "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real loggammapitau;
protected
  Real log_apitau = calc_log_apitau(
                              T,p,X);
algorithm
  loggammapitau :=log_apitau;
annotation(smoothOrder=20);
end calc_loggammapitau;
