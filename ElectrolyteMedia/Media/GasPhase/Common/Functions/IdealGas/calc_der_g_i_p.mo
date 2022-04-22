within ElectrolyteMedia.Media.GasPhase.Common.Functions.IdealGas;
function calc_der_g_i_p
  "Calculates pressure derivative of partial gibbs free energy of gas species"
  input SI.Temperature T "Temperature";
  input SI.Pressure p "Pressure";
  output Real[nGfun] gdp(unit="J/(kg.Pa)") "Specific Gibbs free energy at temperature T";
algorithm
  for i in 1:nGfun loop
    gdp[i] := datafun[i].R*T/p;
  end for;
  annotation (Inline=false,smoothOrder=2);
end calc_der_g_i_p;
