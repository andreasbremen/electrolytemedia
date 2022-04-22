within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson;
function calc_beta
  "Isobaric expansion coefficient with Peng Robinson EOS"
  input SI.Temperature T;
  input SI.Density d;
  input SI.MassFraction[nGfun] X;

  output Real beta(unit="1/K");
protected
  Real v_dT(unit="m3/(kg.K)") = calc_der_v_T(
    T,
    d,
    X);
  SI.SpecificVolume v=calc_v(d);

algorithm
  beta :=1/v*v_dT;
end calc_beta;
