within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson.Molar;
function calc_z
  "Calculates compressibility factor of gas mixture with Peng Robinson EOS"

  input SI.Temperature T;
  input SI.Density d;
  input SI.MoleFraction[nGfun] y_i;

  output Real z(unit="");
protected
  SI.MolarMass MM = Molar.calc_MM(y_i);
  SI.MolarVolume v = MM/d;
  SI.Pressure p=Molar.calc_p(T,d,y_i);
algorithm
  z := p*v/(Modelica.Constants.R*T);
end calc_z;
