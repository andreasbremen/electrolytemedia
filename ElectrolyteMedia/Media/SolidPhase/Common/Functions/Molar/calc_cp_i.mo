within ElectrolyteMedia.Media.SolidPhase.Common.Functions.Molar;
function calc_cp_i "Calculates specific molar heat capacity of mineral species"

  input SI.Temperature T;

  output SI.MolarHeatCapacity c_p[nSfun];
algorithm
  c_p := datafun[:].a + datafun[:].b*T + datafun[:].c*T^(-2) + datafun[:].d*T^(-0.5);
end calc_cp_i;
