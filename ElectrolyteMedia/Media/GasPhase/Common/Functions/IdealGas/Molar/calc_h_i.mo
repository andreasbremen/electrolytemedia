within ElectrolyteMedia.Media.GasPhase.Common.Functions.IdealGas.Molar;
function calc_h_i "Molar enthalpy at temperature T"
  input SI.Temperature T "Temperature";
  output SI.MolarEnthalpy[nGfun] h "Molar enthalpy at temperature T";
protected
  SI.MolarEnthalpy[nGfun] int_cp=calc_int_cp_i_T(T);
algorithm
  for i in 1:nGfun loop
    h[i] := datafun[i].H_ref + int_cp[i];
  end for;
  annotation (Inline=false,smoothOrder=2, derivative=calc_der_h_i_T);
end calc_h_i;
