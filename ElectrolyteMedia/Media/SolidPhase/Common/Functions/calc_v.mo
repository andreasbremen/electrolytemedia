within ElectrolyteMedia.Media.SolidPhase.Common.Functions;
function calc_v "Calculates specific volume of a mixture of mineral species"

  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nSfun] X;
  output SI.SpecificVolume v_mix;

protected
  SI.SpecificVolume[nSfun] v=calc_v_i(T, p);

algorithm

  v_mix :=v*X;

  annotation(smoothOrder=5);
end calc_v;
