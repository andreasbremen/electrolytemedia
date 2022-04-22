within ElectrolyteMedia.Media.GasPhase.Common.Functions.IdealGas;
function calc_der_g_i_T
  "Calculates temperature derivative of partial gibbs free energy of gas species"
  input SI.Temperature T "Temperature";
  output SI.SpecificEntropy[nGfun] gdT "Specific Gibbs free energy at temperature T";
algorithm
  gdT :=-calc_s_i(T);
  annotation (Inline=false,smoothOrder=2);
end calc_der_g_i_T;
