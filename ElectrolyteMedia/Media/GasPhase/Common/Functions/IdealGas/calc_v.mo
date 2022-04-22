within ElectrolyteMedia.Media.GasPhase.Common.Functions.IdealGas;
function calc_v "Calculates ideal gas specific volume"
  input SI.Temperature T "Temperature";
  input SI.Pressure p "Pressure";
  input SI.MassFraction[nGfun] X "Mass fraction";
  output SI.SpecificVolume v "Specific enthalpy at temperature T";
algorithm
  v := (sum(X[i]*datafun[i].R for i in 1:nGfun)*T)/p;
end calc_v;
