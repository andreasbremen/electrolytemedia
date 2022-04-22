within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson.Molar;
function calc_h
  "Calculates molar enthalpy of gas phase with Peng Robinson"
  input SI.Temperature T;
  input SI.Density d;
  input SI.MoleFraction[nGfun] y_i;

  output SI.MolarEnthalpy h;
protected
  SI.MolarMass MM = Molar.calc_MM(y_i);
  SI.MolarVolume v = MM/d;
  SI.Pressure p = Molar.calc_p(T,d,y_i);
  SI.MolarEnergy u=Molar.calc_u(T,d,y_i);
algorithm
  h := u+p*v;

end calc_h;
