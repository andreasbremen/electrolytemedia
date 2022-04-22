within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solute.tautau;
function calc_f_gammatautau "Helper function to calculate Gibbs derivative"

  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real f_gammatautau;
protected
  parameter Real b = 1.2;
  Real A_phitautau=DebyeHueckel.Reduced.Additional.d2g_dtau2.calc_d2Aln_dtau2(T,p);
  Real I = calc_I(X);
algorithm
  f_gammatautau := -A_phitautau * ( I ^ 0.5 / ( 1 + b * I ^ 0.5)  + ( 2 / b)  * log(  1 + b * I ^ 0.5));
end calc_f_gammatautau;
