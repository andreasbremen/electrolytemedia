within ElectrolyteMedia.Media.SolidPhase.Common.Functions.Molar;
function calc_h_i "Calculates enthalpy of mineral species"

  input SI.Temperature T;
  input SI.Pressure p;

  output SI.MolarEnthalpy h[nSfun];
protected
  SI.MolarEnergy vdp[nSfun]=calc_int_v_p(T, p);
  SI.MolarEnergy int_c_p_T[nSfun]=calc_int_cp_T(T);
algorithm
  h := datafun[:].H_ref + int_c_p_T + vdp;
  annotation(smoothOrder=5);
end calc_h_i;
