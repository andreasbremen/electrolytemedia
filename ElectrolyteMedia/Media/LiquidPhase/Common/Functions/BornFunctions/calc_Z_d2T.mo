within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.BornFunctions;
function calc_Z_d2T "Calculates Born function X"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  output Real X;
protected
  Real eps=calc_eps(T, p);
  Real eps_d2T=calc_der_eps_TT(T, p);
  Real Y=calc_Z_dT(T, p);
algorithm
  X := 1 / eps ^ 2 * eps_d2T - 2 * eps * Y ^ 2;
  annotation(smoothOrder=5);
end calc_Z_d2T;
