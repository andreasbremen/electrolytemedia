within ElectrolyteMedia.Media.SolidPhase.Common.Functions.Molar;
function calc_u_i
  "Calculates molar internal energy of mineral species"

  input SI.Temperature T;
  input SI.Pressure p;

  output SI.MolarEnthalpy u_i[nSfun];
protected
  SI.MolarVolume v_i[nSfun]=calc_v_i(T, p);
  SI.MolarEnthalpy h_i[nSfun]=calc_h_i(T, p);
algorithm
  u_i :=h_i - p*v_i;
  annotation(smoothOrder=5);
end calc_u_i;
