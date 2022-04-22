within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.BornFunctions;
function calc_Z_dT "Calculates Born function Y"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  output Real Y;
protected
  Real eps=calc_eps(T, p);
  Real eps_dT=calc_der_eps_T(T, p);
algorithm
  Y := 1 / eps ^ 2 * eps_dT;
  annotation(smoothOrder=5);
end calc_Z_dT;
