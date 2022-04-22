within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.Solutes.Reduced.Sym;
function calc_tau "Calculates dimensionless temperature for asymmetrical mixing"
  input SI.Temperature T;
  output Real tau;
algorithm
  tau :=Tref/T;
end calc_tau;
