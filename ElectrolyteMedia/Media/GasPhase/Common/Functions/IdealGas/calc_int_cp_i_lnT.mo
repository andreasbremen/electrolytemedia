within ElectrolyteMedia.Media.GasPhase.Common.Functions.IdealGas;
function calc_int_cp_i_lnT "Specific heat capacity integral over lnT"
  input SI.Temperature T "Temperature";
  output SI.SpecificEntropy[nGfun] int_cp_lnT "Specific heat capacity integral over lnT";
algorithm
  int_cp_lnT :=Molar.calc_int_cp_i_lnT(T) ./ datafun.MM;
end calc_int_cp_i_lnT;
