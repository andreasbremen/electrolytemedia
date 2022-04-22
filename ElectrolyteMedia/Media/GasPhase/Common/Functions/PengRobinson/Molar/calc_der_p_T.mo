within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson.Molar;
function calc_der_p_T
  "Calculates derivative w.r.t. T at const. d of pressure of gas mixture with Peng Robinson EOS"

  input SI.Temperature T;
  input SI.Density d;
  input SI.MoleFraction[nGfun] y;

  output Real p_dT(unit="Pa/K");
protected
  SI.MolarMass MM = Molar.calc_MM(y);
  SI.MolarVolume v = MM/d;
  SI.MolarVolume b=Molar.calc_b(y);
  Real a_alpha_dT(unit="N.m4/(mol2.K)") = Molar.calc_a_dT(T, y);
algorithm
  p_dT :=Modelica.Constants.R/(v - b) - a_alpha_dT/(v*(v+b) +b*(v-b));
end calc_der_p_T;
