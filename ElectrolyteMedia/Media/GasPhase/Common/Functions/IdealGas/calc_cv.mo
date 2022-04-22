within ElectrolyteMedia.Media.GasPhase.Common.Functions.IdealGas;
function calc_cv
  "Calculates specific heat capacity at constant volume of gas phase with ideal gas model"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nGfun] X;

  output SI.SpecificHeatCapacityAtConstantVolume cv;
protected
  SI.SpecificHeatCapacityAtConstantPressure[nGfun] cpi=calc_cp_i(T);
algorithm
  cv :=X*cpi - sum(datafun[i].R*X[i] for i in 1:nGfun);

end calc_cv;
