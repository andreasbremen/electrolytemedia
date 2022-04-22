within ElectrolyteMedia.Media.SolidPhase.Common.Functions.Molar;
function calc_s_i "Calculates thermodynamic properties of mineral species"

  input SI.Temperature T;

  output SI.MolarEntropy s[nSfun];
protected
  Real int_c_p_lnT[nSfun];
algorithm
  int_c_p_lnT :=calc_int_cp_lnT(T);
  s := datafun[:].S_ref + int_c_p_lnT;
end calc_s_i;
