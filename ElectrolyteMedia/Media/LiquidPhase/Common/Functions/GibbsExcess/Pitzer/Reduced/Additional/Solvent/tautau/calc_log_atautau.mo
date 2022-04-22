within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solvent.tautau;
function calc_log_atautau "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real log_atautau;
protected
  Real ln_atautau = calc_ln_atautau(T,p,X);
algorithm
  log_atautau :=ln_atautau*log10(exp(1));

end calc_log_atautau;
