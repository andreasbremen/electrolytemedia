within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solvent.pi;
function calc_loggammapi "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real loggammapi;
protected
  Real log_api = calc_log_api(T,p,X);
algorithm
  loggammapi :=log_api;
annotation(smoothOrder=20);
end calc_loggammapi;
