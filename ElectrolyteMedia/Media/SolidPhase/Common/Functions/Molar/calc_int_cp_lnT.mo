within ElectrolyteMedia.Media.SolidPhase.Common.Functions.Molar;
function calc_int_cp_lnT
  "Calculates integral of specific molar heat capacity of mineral species over ln(T)"

  input SI.Temperature T;

  output SI.MolarEntropy int_c_p_lnT[nSfun];
algorithm
  int_c_p_lnT := datafun[:].a*log(T/T_0) + datafun[:].b*(T - T_0) - datafun[:].c/2
    *(1/T^2 - 1/T_0^2) - datafun[:].d/0.5*(1/T^0.5 - 1/T_0^0.5);
  annotation (smoothOrder=20);
end calc_int_cp_lnT;
