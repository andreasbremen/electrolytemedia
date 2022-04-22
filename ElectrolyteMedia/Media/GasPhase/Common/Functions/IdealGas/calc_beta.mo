within ElectrolyteMedia.Media.GasPhase.Common.Functions.IdealGas;
function calc_beta "Isobaric expansion coefficient"
  input SI.Temperature T;
  output Real beta(unit="1/K");
algorithm
  beta :=1/T;
end calc_beta;
