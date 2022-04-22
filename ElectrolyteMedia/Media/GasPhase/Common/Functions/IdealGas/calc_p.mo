within ElectrolyteMedia.Media.GasPhase.Common.Functions.IdealGas;
function calc_p "Calculates ideal gas pressure"
  input SI.Temperature T "Temperature";
  input SI.Density d "Density";
  input SI.MassFraction[nGfun] X "Mass fraction";
  output SI.Pressure p "Pressure";
protected
  SI.SpecificVolume v = 1/d;
algorithm
  p := (sum(X[i]*datafun[i].R for i in 1:nGfun)*T)/v;
end calc_p;
