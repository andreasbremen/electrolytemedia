within ElectrolyteMedia.Media.GasPhase.Common.Functions.IdealGas.Molar;
function calc_der_h_i_T "Derivative of molar enthalpy at temperature T"
  input SI.Temperature T "Temperature";
  input Real dT( unit="K.s-1");
  output Real[nGfun] h_der(unit="J.mol-1.s-1") "Derivative of molar enthalpy at temperature T";
protected
  SI.MolarHeatCapacity[nGfun] cp=calc_cp_i(T);
algorithm
  h_der := dT*cp;
  annotation (Inline=false,smoothOrder=2);
end calc_der_h_i_T;
