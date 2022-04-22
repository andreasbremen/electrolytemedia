within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.Solutes.Reduced;
function calc_pi "Calculates dimensionless pressure"
  input SI.Pressure p;
  output Real pi;
algorithm
  pi :=p/pref;
end calc_pi;
