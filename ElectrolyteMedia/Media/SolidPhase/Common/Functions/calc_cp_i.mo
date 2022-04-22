within ElectrolyteMedia.Media.SolidPhase.Common.Functions;
function calc_cp_i
  "Calculates specific molar heat capacity of mineral species"

  input SI.Temperature T;

  output SI.SpecificHeatCapacity c_p[nSfun];

protected
  SI.MolarHeatCapacity c_p_molar[nSfun]=Molar.calc_cp_i(T);
algorithm

  for i in 1:nSfun loop
    c_p[i] := c_p_molar[i]/datafun[i].MM;
  end for;

end calc_cp_i;
