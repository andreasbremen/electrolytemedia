within ElectrolyteMedia.Media.SolidPhase.Common.Functions;
function calc_XfromPhi
  "calculate mass fraction vector from volume fraction vector"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.VolumeFraction[nSfun] Phi;
  output SI.MassFraction[nSfun] X;
protected
  SI.SpecificVolume[nSfun] v_i = calc_v_i(T,p);
  SI.Density[nSfun] d_i = Phi./v_i;
algorithm
  X :=d_i/sum(d_i);

end calc_XfromPhi;
