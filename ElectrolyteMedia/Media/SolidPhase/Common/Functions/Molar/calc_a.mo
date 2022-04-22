within ElectrolyteMedia.Media.SolidPhase.Common.Functions.Molar;
function calc_a "Calculates a of mineral EOS from Hooland & Powell 2011"

  output Real a[nSfun];
algorithm
  a := (ones(nSfun) + datafun[:].k_0_) ./ (ones(nSfun) + datafun[:].k_0_ + datafun[:].k_0 .* datafun[:].k_0__);
  annotation(smoothOrder=5);
end calc_a;
