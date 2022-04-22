within ElectrolyteMedia.Media.GasPhase.Common.Functions.IdealGas.Molar;
function calc_s0_i "Molar entropy at temperature T"
  input SI.Temperature T "Temperature";
  output SI.MolarEntropy[nGfun] s;

protected
  SI.MolarEntropy[nGfun] int_cp=calc_int_cp_i_lnT(T);
algorithm
  for i in 1:nGfun loop
    s[i] := datafun[i].S_ref + int_cp[i];
  end for;
end calc_s0_i;
