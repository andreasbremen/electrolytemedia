within ElectrolyteMedia.Media.GasPhase.Common.Functions.IdealGas;
function calc_kappa "Isothermal compressibility"
  input SI.Pressure p;
  output Real kappa(unit="1/Pa");
algorithm
  kappa :=1/p;
end calc_kappa;
