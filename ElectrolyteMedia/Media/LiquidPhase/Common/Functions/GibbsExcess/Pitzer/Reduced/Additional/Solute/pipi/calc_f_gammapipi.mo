within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solute.pipi;
function calc_f_gammapipi "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real f_gammapipi;
protected
  parameter Real b = 1.2;
  Real A_phipipi=GibbsExcess.DebyeHueckel.Reduced.Additional.d2g_dpi2.calc_d2Aln_dpi2(T,p);
  Real I = calc_I(X);
algorithm
  f_gammapipi := -A_phipipi * ( I ^ 0.5 / ( 1 + b * I ^ 0.5)  + ( 2 / b)  * log(  1 + b * I ^ 0.5));
end calc_f_gammapipi;
