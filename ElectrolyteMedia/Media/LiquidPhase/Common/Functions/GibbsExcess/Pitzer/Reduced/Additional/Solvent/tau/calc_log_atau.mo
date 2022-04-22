within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solvent.tau;
function calc_log_atau "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real log_a;
protected
  Real ln_atau = calc_ln_atau(T,p,X);
algorithm
  log_a :=ln_atau*log10(exp(1));

end calc_log_atau;
