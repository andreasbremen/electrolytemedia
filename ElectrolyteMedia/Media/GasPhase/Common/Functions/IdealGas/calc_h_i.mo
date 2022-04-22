within ElectrolyteMedia.Media.GasPhase.Common.Functions.IdealGas;
function calc_h_i "Specific enthalpy at temperature T"
  input SI.Temperature T "Temperature";
  output SI.SpecificEnthalpy[nGfun] h "Specific enthalpy at temperature T";
algorithm
  h :=Molar.calc_h_i(T) ./ datafun.MM;
  annotation (Inline=false,smoothOrder=2);
end calc_h_i;
