within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solute.pi;
function calc_f_gammapi "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real f_gammapi;
protected
  parameter Real b = 1.2;
  Real A_phipi=GibbsExcess.DebyeHueckel.Reduced.Additional.dg_dpi.calc_dAln_dpi(T,p);
  Real I = calc_I(X);
algorithm
  f_gammapi := -A_phipi * ( I ^ 0.5 / ( 1 + b * I ^ 0.5)  + ( 2 / b)  * log(  1 + b * I ^ 0.5));
end calc_f_gammapi;
