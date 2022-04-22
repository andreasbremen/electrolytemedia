within ElectrolyteMedia.Media.SolidPhase.Common.Functions;
function calc_s_i
  "Calculates specific entropy of mineral species"

  input SI.Temperature T;

  output SI.SpecificEntropy s[nSfun];

protected
  SI.MolarEntropy[nSfun] s_m=Molar.calc_s_i(T);

algorithm

  for i in 1:nSfun loop
    s[i] :=s_m[i]/datafun[i].MM;
  end for;

end calc_s_i;
