within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson;
function calc_der_p_T
  "Calculates derivative w.r.t. T at const. d of pressure of gas mixture with Peng Robinson EOS"

  input SI.Temperature T;
  input SI.Density d;
  input SI.MassFraction[nGfun] X;

  output Real p_dT(unit="Pa/K");
protected
  SI.MoleFraction[nGfun] y=calc_Y(X);
algorithm
  p_dT :=Molar.calc_der_p_T(
      T,
      d,
      y);
end calc_der_p_T;
