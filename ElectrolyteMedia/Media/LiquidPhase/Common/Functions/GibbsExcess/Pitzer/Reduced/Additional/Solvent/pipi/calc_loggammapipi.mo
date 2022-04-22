within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solvent.pipi;
function calc_loggammapipi "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real loggammapipi;
protected
  Real log_apipi = calc_log_apipi(T,p,X);
algorithm
  loggammapipi :=log_apipi;
annotation(smoothOrder=20);
end calc_loggammapipi;
