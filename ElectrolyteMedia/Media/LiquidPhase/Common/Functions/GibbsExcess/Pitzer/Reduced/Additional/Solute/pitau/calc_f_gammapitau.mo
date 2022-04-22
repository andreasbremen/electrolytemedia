within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solute.pitau;
function calc_f_gammapitau "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real f_gammapitau;
protected
  parameter Real b = 1.2;
  Real A_phipitau=GibbsExcess.DebyeHueckel.Reduced.Additional.d2g_dtaudpi.calc_d2Aln_dtaudpi(T,p);
  Real I = calc_I(X);
algorithm
  f_gammapitau := -A_phipitau * ( I ^ 0.5 / ( 1 + b * I ^ 0.5)  + ( 2 / b)  * log(  1 + b * I ^ 0.5));
end calc_f_gammapitau;
