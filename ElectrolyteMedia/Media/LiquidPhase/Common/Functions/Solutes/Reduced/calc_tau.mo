within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.Solutes.Reduced;
function calc_tau "Calculates dimensionless temperature"
  input SI.Temperature T;
  output Real tau;
algorithm
  tau :=Tref/T;
end calc_tau;
