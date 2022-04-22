within ElectrolyteMedia.Media.SolidPhase.Common.Functions;
function calc_u "Calculates specific internal energy of a mixture of mineral species"

  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nSfun] X;
  output SI.SpecificEnthalpy u_mix;

protected
  SI.SpecificEnthalpy u_i[nSfun] = calc_u_i(T,p);

algorithm

  u_mix :=u_i*X;

  annotation(smoothOrder=5);
end calc_u;
