within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solvent.pitau;
function calc_log_apitau "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real log_apitau;
protected
  Real ln_apitau = calc_ln_apitau(
                            T,p,X);
algorithm
  log_apitau :=ln_apitau*log10(exp(1));

end calc_log_apitau;
