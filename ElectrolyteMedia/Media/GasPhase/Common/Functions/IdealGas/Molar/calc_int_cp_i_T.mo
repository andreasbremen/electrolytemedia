within ElectrolyteMedia.Media.GasPhase.Common.Functions.IdealGas.Molar;
function calc_int_cp_i_T "Molar heat capacity integral over T"
  input SI.Temperature T "Temperature";
  output SI.MolarEnergy[nGfun] int_cp_T "Molar heat capacity integral over T";
algorithm
  for i in 1:nGfun loop
    int_cp_T[i] := datafun[i].a*(T - T0) + 0.5*datafun[i].b*(T^2 - T0^2) - datafun[i].c*(T^(-1)- T0^(-1)) + 2*datafun[i].d*(T^0.5 - T0^0.5);
  end for;
  annotation (smoothOrder=5);
end calc_int_cp_i_T;
