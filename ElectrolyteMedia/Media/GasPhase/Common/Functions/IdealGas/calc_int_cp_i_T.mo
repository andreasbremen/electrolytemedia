within ElectrolyteMedia.Media.GasPhase.Common.Functions.IdealGas;
function calc_int_cp_i_T "Specific heat capacity integral over T"
  input SI.Temperature T "Temperature";
  output SI.SpecificEnergy[nGfun] int_cp_T "Specific heat capacity integral over T";
algorithm
  int_cp_T :=Molar.calc_int_cp_i_T(T) ./ datafun.MM;
end calc_int_cp_i_T;
