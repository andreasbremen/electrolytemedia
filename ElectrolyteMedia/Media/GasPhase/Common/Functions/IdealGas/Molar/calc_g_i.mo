within ElectrolyteMedia.Media.GasPhase.Common.Functions.IdealGas.Molar;
function calc_g_i "Molar Gibbs free energy at temperature T"
  input SI.Temperature T "Temperature";
  output SI.MolarEnergy[nGfun] g "Molar Gibbs free energy at temperature T";
protected
  SI.MolarEnergy[nGfun] intcpT=calc_int_cp_i_T(T);
  SI.MolarEntropy[nGfun] intcplnT=calc_int_cp_i_lnT(T);
algorithm
  for i in 1:nGfun loop
    g[i] :=datafun[i].G_ref + intcpT[i] - T*intcplnT[i] - datafun[i].S_ref*(T - T0);
  end for;
  annotation(smoothOrder=5);
end calc_g_i;
