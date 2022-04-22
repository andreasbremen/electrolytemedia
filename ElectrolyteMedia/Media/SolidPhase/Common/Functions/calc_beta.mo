within ElectrolyteMedia.Media.SolidPhase.Common.Functions;
function calc_beta
  "Isobaric expansion coefficient"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nSfun] X;
  output Real beta(unit="1/K");

protected
  Real v_dT=calc_der_v_T(
      T,
      p,
      X);
  SI.SpecificVolume v=calc_v(T,p,X);

algorithm
  beta :=1/v*v_dT;

end calc_beta;
