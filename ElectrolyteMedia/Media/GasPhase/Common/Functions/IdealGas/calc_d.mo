within ElectrolyteMedia.Media.GasPhase.Common.Functions.IdealGas;
function calc_d "Calculates ideal gas density"
  input SI.Temperature T "Temperature";
  input SI.Pressure p "Pressure";
  input SI.MassFraction[nGfun] X "Mass fraction";
  output SI.Density d "Density at temperature T and pressure p";
protected
  SI.SpecificVolume v = calc_v(T,p,X);
algorithm
  d := 1/v;
end calc_d;
