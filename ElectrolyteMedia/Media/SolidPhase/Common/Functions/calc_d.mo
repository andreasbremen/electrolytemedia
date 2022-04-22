within ElectrolyteMedia.Media.SolidPhase.Common.Functions;
function calc_d "Calculates density of a mixture of mineral species"

  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nSfun] X;
  output SI.Density d_mix;

algorithm

  d_mix :=1/calc_v(
      T,
      p,
      X);

end calc_d;
