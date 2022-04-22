within ElectrolyteMedia.Media.GasPhase.Common.Functions.IdealGas.Molar;
function calc_int_cp_i_lnT "Molar heat capacity integral over lnT"
  input SI.Temperature T "Temperature";
  output SI.MolarEntropy[nGfun] int_cp_lnT "Molar heat capacity integral over lnT";
algorithm
  for i in 1:nGfun loop
    int_cp_lnT[i] := datafun[i].a * log(T / T0) + datafun[i].b * (T - T0) - datafun[i].c / 2 * (1 / T ^ 2 - 1 / T0 ^ 2) - datafun[i].d / 0.5 * (1 / T ^ 0.5 - 1 / T0 ^ 0.5);
  end for;
  annotation(smoothOrder=5);
end calc_int_cp_i_lnT;
