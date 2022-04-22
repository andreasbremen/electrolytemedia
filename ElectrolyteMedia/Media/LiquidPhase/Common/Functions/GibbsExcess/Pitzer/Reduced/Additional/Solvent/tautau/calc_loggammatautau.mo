within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solvent.tautau;
function calc_loggammatautau "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real loggammatautau;
protected
  Real log_atautau = calc_log_atautau(T,p,X);
algorithm
  loggammatautau :=log_atautau;
annotation(smoothOrder=20);
end calc_loggammatautau;
