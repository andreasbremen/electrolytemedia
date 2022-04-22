within ElectrolyteMedia.Media.LiquidPhase.Common.MixtureLiquid;
function calc_mred_from_mfull
  "Calculates reduced mass vector from full mass vector"

  input Real[nF] mfull;
  output Real[nX] mred;

algorithm
  mred :=transpose(lambda_mass)*mfull;

end calc_mred_from_mfull;
