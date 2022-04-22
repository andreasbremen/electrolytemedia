within ElectrolyteMedia.Media.SolidPhase.Common.Functions;
function calc_der_d_p
  "Calculates pressure derivative of density of solid mixture at T and p"

  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nSfun] X;

  output SI.DerDensityByPressure ddpT;
algorithm
  ddpT :=0;

end calc_der_d_p;
