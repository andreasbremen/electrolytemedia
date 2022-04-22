within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson.Molar;
function calc_B_dT_p
  "Calculates derivative B w.r.t. T at const p in Peng Robinson EOS"

  input SI.Temperature T;
  input SI.Density d;
  input SI.MoleFraction[nGfun] y_i;

  output Real B_dT(unit="1/K");
protected
  SI.Pressure p=Molar.calc_p(T,d,y_i);
  SI.MolarVolume b=Molar.calc_b(y_i);
algorithm
  B_dT := -b*p/(Modelica.Constants.R*T^2);
end calc_B_dT_p;
