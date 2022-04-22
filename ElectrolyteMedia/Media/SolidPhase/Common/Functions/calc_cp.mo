within ElectrolyteMedia.Media.SolidPhase.Common.Functions;
function calc_cp "Calculates specific molar heat capacity of a mixture of mineral species"

  input SI.Temperature T;
  input SI.MassFraction[nSfun] X;
  output SI.SpecificHeatCapacityAtConstantPressure cp_mix;

protected
  SI.SpecificHeatCapacityAtConstantPressure cp[nSfun]=calc_cp_i(T);

algorithm

  cp_mix :=cp*X;

end calc_cp;
