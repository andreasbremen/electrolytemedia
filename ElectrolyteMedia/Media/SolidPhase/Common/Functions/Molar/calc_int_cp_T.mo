within ElectrolyteMedia.Media.SolidPhase.Common.Functions.Molar;
function calc_int_cp_T
  "Calculates integral of specific molar heat capacity of mineral species over T"

  input SI.Temperature T;

  output SI.MolarEnergy int_c_p_dT[nSfun];
algorithm
  int_c_p_dT := datafun[:].a*(T - T_0) + 0.5*datafun[:].b*(T^2 - T_0^2) - datafun[:].c*(
    T^(-1) - T_0^(-1)) + 2*datafun[:].d*(T^0.5 - T_0^0.5);
  annotation (smoothOrder=20);
end calc_int_cp_T;
