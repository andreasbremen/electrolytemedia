within ElectrolyteMedia.Media.GasPhase.Common.Functions.IdealGas.Molar;
function calc_cp_i
  "Compute molar heat capacity at constant pressure from temperature and gas data"
  input SI.Temperature T "Temperature";
  output SI.MolarHeatCapacity[nGfun] cp "Specific heat capacity at temperature T";
algorithm
  for i in 1:nGfun loop
    cp[i] := datafun[i].a + datafun[i].b * T + datafun[i].c * T ^ (-2) + datafun[i].d * T ^ (-0.5);
  end for;
  annotation (smoothOrder=5);
end calc_cp_i;
