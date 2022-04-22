within ElectrolyteMedia.Media.GasPhase.Common.Functions.IdealGas;
function calc_s_i "Specific entropy at temperature T"
  input SI.Temperature T "Temperature";
  output SI.SpecificEntropy[nGfun] s;
algorithm
  s :=Molar.calc_s0_i(T) ./ datafun.MM;
end calc_s_i;
