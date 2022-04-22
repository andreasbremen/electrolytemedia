within ElectrolyteMedia.Media.GasPhase.Common.Functions.IdealGas;
function calc_cp_i
  "Compute specific heat capacity at constant pressure from temperature and gas data"
  input SI.Temperature T "Temperature";
  output SI.SpecificHeatCapacity[nGfun] cp "Specific heat capacity at temperature T";
algorithm
  cp :=Molar.calc_cp_i(T) ./ datafun.MM;
end calc_cp_i;
