within ElectrolyteMedia.Media.SolidPhase.Common.Functions.Molar;
function calc_g_i "Calculates Gibbs free energy of mineral species"

  input SI.Temperature T;
  input SI.Pressure p;

  output SI.MolarEnergy g[nSfun];
protected
  SI.MolarEnergy vdp[nSfun]=calc_int_v_p(T, p);
  SI.MolarEntropy int_c_p_lnT[nSfun]=calc_int_cp_lnT(T);
  SI.MolarEnergy int_c_p_dT[nSfun]=calc_int_cp_T(T);
algorithm
  g := datafun[:].G_ref + vdp + int_c_p_dT - T * int_c_p_lnT - datafun[:].S_ref * (T - T_0);
  annotation(smoothOrder=5);
end calc_g_i;
