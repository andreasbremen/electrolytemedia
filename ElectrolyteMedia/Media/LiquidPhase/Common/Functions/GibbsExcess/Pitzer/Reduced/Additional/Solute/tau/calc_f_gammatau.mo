within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solute.tau;
function calc_f_gammatau "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real f_gammatau;
protected
  parameter Real b = 1.2;
  Real A_phitau=DebyeHueckel.Reduced.Additional.dg_dtau.calc_dAln_dtau(T,p);
  Real I = calc_I(X);
algorithm
  f_gammatau := -A_phitau * ( I ^ 0.5 / ( 1 + b * I ^ 0.5)  + ( 2 / b)  * log(  1 + b * I ^ 0.5));
end calc_f_gammatau;
