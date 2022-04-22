within ElectrolyteMedia.Media.GasPhase.Common.Functions.IdealGas;
function calc_der_h_i_T "Derivative of specific enthalpy at temperature T"
  input SI.Temperature T "Temperature";
  input Real dT( unit="K.s-1");
  output Real[nGfun] h_der(unit="J.kg-1.s-1") "Derivative of specific enthalpy at temperature T";
algorithm
  h_der :=Molar.calc_der_h_i_T(T, dT) ./ datafun.MM;
  annotation (Inline=false,smoothOrder=2);
end calc_der_h_i_T;
