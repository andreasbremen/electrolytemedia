within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.BornFunctions;
function calc_Z_dp "Calculates Born function Q"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  output Real Q;
protected
  Real eps=calc_eps(T, p);
  Real eps_dp=calc_der_eps_p(T, p);
algorithm
  Q := 1 / eps ^ 2 * eps_dp;
  annotation(smoothOrder=5);
end calc_Z_dp;
