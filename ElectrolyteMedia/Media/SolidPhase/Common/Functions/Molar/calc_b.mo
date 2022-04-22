within ElectrolyteMedia.Media.SolidPhase.Common.Functions.Molar;
function calc_b "Calculates b of mineral EOS from Hooland & Powell 2011"

  output Real b[nSfun];
algorithm
  b := datafun[:].k_0_ ./ datafun[:].k_0 - datafun[:].k_0__ ./ (ones(nSfun) + datafun[:].k_0_);
  annotation(smoothOrder=5);
end calc_b;
