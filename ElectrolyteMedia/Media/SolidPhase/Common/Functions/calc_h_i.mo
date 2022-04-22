within ElectrolyteMedia.Media.SolidPhase.Common.Functions;
function calc_h_i   "Calculates specific enthalpy of mineral species"

  input SI.Temperature T;
  input SI.Pressure p;

  output SI.SpecificEnthalpy h[nSfun];

protected
  SI.MolarEnthalpy[nSfun] h_m=Molar.calc_h_i(T, p);

algorithm

  for i in 1:nSfun loop
    h[i] :=h_m[i]/datafun[i].MM;
  end for;

end calc_h_i;
