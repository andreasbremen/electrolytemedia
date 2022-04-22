within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.BornFunctions;
function calc_Z "Calculates Born function Z"
  input SI.Temperature T;
  input SI.Pressure p;
  output Real Z;
protected
  Real eps=calc_eps(T, p);
algorithm
  Z := -1 / eps;
  annotation(smoothOrder=5);
end calc_Z;
