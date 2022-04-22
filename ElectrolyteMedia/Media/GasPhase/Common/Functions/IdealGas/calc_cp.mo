within ElectrolyteMedia.Media.GasPhase.Common.Functions.IdealGas;
function calc_cp
  "Calculates specific heat capacity at constant pressure of gas phase with ideal gas model"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nGfun] X;

  output SI.SpecificHeatCapacityAtConstantPressure cp;
protected
  SI.SpecificHeatCapacityAtConstantPressure[nGfun] cpi=calc_cp_i(T);
algorithm
  cp :=X*cpi;

end calc_cp;
