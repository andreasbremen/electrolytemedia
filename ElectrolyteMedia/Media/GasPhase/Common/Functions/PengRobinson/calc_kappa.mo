within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson;
function calc_kappa "Isothermal compressibility with Peng Robinson EOS"
  input SI.Temperature T;
  input SI.Density d;
  input SI.MoleFraction[nGfun] X;

  output Real alpha(unit="1/Pa");
protected
  Real v_dp(unit="m3/(kg.Pa)") = calc_der_v_p(
    T,
    d,
    X);
  SI.SpecificVolume v=calc_v(d);

algorithm
  alpha :=1/v*v_dp;
end calc_kappa;
