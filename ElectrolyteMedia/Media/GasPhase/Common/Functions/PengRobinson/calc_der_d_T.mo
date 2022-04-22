within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson;
function calc_der_d_T
  "Calculates temperature derivative of at T and p with Peng-Robinson EOS"
    input SI.Temperature T;
  input SI.Density d;
  input SI.MassFraction[nGfun] X;

  output SI.DerDensityByTemperature ddTp;
protected
  Real vdT(unit="m3/(kg.K)") = calc_der_v_T(
    T,
    d,
    X);
algorithm
  ddTp := -d*d*vdT;

end calc_der_d_T;
