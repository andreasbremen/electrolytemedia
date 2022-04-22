within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson;
function calc_der_d_p
  "Calculates pressure derivative of density at T and p with Peng-Robinson EOS"
    input SI.Temperature T;
  input SI.Density d;
  input SI.MassFraction[nGfun] X;

  output SI.DerDensityByPressure ddpT;
protected
  Real vdp(unit="m3/(kg.Pa)") = calc_der_v_p(
    T,
    d,
    X);
algorithm
  ddpT := -d*d*vdp;

end calc_der_d_p;
