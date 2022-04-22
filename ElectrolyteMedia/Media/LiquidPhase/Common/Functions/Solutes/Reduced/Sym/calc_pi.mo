within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.Solutes.Reduced.Sym;
function calc_pi "Calculates dimensionless pressure for asymmetrical mixing"
  input SI.Pressure p;
  output Real pi;
algorithm
  pi :=p/pref;
end calc_pi;
