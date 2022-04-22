within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson;
function calc_der_v_p
  "Calculates derivative w.r.t. p at const. T of specific volume of gas mixture with Peng Robinson EOS"
    input SI.Temperature T;
  input SI.Density d;
  input SI.MassFraction[nGfun] X;

  output Real v_dp(unit="m3/(kg.Pa)");
protected
  Real v_dT(unit="m3/(kg.K)") = calc_der_v_T(
    T,
    d,
    X);
  Real p_dT(unit="Pa/K") = calc_der_p_T(
    T,
    d,
    X);
algorithm
  v_dp := -v_dT/p_dT;

end calc_der_v_p;
