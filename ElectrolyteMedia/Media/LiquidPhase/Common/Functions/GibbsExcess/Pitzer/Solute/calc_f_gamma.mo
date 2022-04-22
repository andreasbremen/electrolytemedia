within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Solute;
function calc_f_gamma
  "Calculates f_gamma that accounts for long range electrostatic effects dependent on ionic strength, temperature and solute properties"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real f_gamma;
protected
  parameter Real b = 1.2;
  Real A_phi=DebyeHueckel.calc_A_ln(T,p);
  Real I = calc_I(X);
algorithm
  f_gamma := -A_phi * ( I ^ 0.5 / ( 1 + b * I ^ 0.5)  + ( 2 / b)  * log(  1 + b * I ^ 0.5));
end calc_f_gamma;
