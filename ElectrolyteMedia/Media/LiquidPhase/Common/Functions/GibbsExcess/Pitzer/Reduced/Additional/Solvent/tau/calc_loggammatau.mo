within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solvent.tau;
function calc_loggammatau "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real loggammatau;
protected
  Real log_atau = calc_log_atau(T,p,X);
algorithm
  loggammatau :=log_atau;
annotation(smoothOrder=20);
end calc_loggammatau;
