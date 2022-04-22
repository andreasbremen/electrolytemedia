within ElectrolyteMedia.Media.SolidPhase.Common.Functions;
function calc_h "Calculates specific enthalpy of a mixture of mineral species"

  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nSfun] X;
  output SI.SpecificEnthalpy h_mix;

protected
  SI.SpecificEnthalpy h[nSfun]=calc_h_i(T, p);

algorithm

  h_mix :=h*X;

end calc_h;
