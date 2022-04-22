within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Debye.Solvent;
function calc_Phi
  "Calculates osmotic coefficient based on Bye limiting law with integration taken from Robinson and Stokes (1959)"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real Phi;
protected
  Real I=calc_I(X);
  Real A=DebyeHueckel.calc_A_log(T, p);
algorithm
  Phi :=1 - A*sqrt(I)/3;
end calc_Phi;
