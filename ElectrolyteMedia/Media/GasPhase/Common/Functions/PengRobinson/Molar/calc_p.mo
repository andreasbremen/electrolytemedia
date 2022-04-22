within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson.Molar;
function calc_p
  "Calculates pressure of gas mixture with Peng Robinson EOS"

  input SI.Temperature T;
  input SI.Density d;
  input SI.MoleFraction[nGfun] y_i;

  output SI.Pressure p;
protected
  SI.MolarMass MM = Molar.calc_MM(y_i);
  SI.MolarVolume v = MM/d;
  SI.MolarVolume b=Molar.calc_b(y_i);
  Real a_alpha(unit="N.m4/mol2") = Molar.calc_a(T, y_i);
algorithm
  p :=Modelica.Constants.R*T/(v - b) - a_alpha/(v*(v+b) +b*(v-b));
end calc_p;
