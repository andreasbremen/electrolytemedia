within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson.Molar;
function calc_s
  "Calculates entropy of gas phase with Peng Robinson"
  input SI.Temperature T;
  input SI.Density d;
  input SI.MoleFraction[nGfun] y_i;

  output SI.MolarEntropy s;

protected
  SI.MolarMass MM = Molar.calc_MM(y_i);
  SI.MolarVolume v = MM/d;
  SI.Pressure p=Molar.calc_p(T,d,y_i);
  SI.MolarEntropy s_ig=IdealGas.Molar.calc_s_i(T, y_i) - Modelica.Constants.R*log(p/
      pref);
  Real a_dT=Molar.calc_a_dT(T,y_i);
  SI.MolarVolume b=Molar.calc_b(y_i);
  Real B = Molar.calc_B(T,d,y_i);
  Real z = Molar.calc_z(T,d,y_i);
algorithm

 s := s_ig + Modelica.Constants.R*log(z-B) + a_dT/(2*sqrt(2)*b)*log((v+(1+sqrt(2))*b)/(v+(1-sqrt(2))*b));

end calc_s;
