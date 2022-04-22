within ElectrolyteMedia.Media.SolidPhase.Common.Functions;
function calc_s "Calculates specific entropy of a mixture of mineral species"

  input SI.Temperature T;
  input SI.MassFraction[nSfun] X;
  output SI.SpecificEntropy s_mix;

protected
  SI.SpecificEntropy s[nSfun]=calc_s_i(T);

algorithm

  s_mix :=s*X;

end calc_s;
