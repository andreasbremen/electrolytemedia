within ElectrolyteMedia.Media.SolidPhase.Common.Functions.Molar;
function calc_c "Calculates c of mineral EOS from Hooland & Powell 2011"

  output Real c[nSfun];
algorithm
  c := (ones(nSfun) + datafun[:].k_0_ + datafun[:].k_0 .* datafun[:].k_0__) ./ (datafun[:].k_0_ .^ 2 + datafun[:].k_0_ - datafun[:].k_0 .* datafun[:].k_0__);
  annotation(smoothOrder=5);
end calc_c;
