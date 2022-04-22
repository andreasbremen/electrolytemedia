within ElectrolyteMedia.Media.SolidPhase.Common.Functions;
function calc_v_i   "Calculates specific volume of mineral species"

  input SI.Temperature T;
  input SI.Pressure p;
  output SI.SpecificVolume v[nSfun];

protected
  SI.MolarVolume[nSfun] v_m=Molar.calc_v_i(T, p);

algorithm

  for i in 1:nSfun loop
    v[i] :=v_m[i]/datafun[i].MM;
  end for;

end calc_v_i;
