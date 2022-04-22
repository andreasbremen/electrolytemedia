within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson.Molar;
function calc_v "Calculates molar volume of a mixture"

  input SI.Density d;
  input SI.MoleFraction y_i[nGfun];
  output SI.MolarVolume v;

protected
  SI.MolarMass MM = Molar.calc_MM(y_i);
algorithm
  v := MM/d;

end calc_v;
