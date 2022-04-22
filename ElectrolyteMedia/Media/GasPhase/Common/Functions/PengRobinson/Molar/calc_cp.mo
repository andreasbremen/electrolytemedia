within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson.Molar;
function calc_cp
  "Calculates molar heat capacity of gas phase with Peng Robinson"
  input SI.Temperature T;
  input SI.Density d;
  input SI.MoleFraction[nGfun] y_i;

  output SI.MolarHeatCapacity cp;
protected
  SI.MolarMass MM = Molar.calc_MM(y_i);
  SI.MolarVolume v = MM/d;
  SI.Pressure p = Molar.calc_p(T,d,y_i);
  SI.MolarHeatCapacity[nGfun] cp_ig_i=IdealGas.Molar.calc_cp_i(T);
  SI.MolarHeatCapacity cp_ig = cp_ig_i*y_i;
  Real z( unit = "") = Molar.calc_z(T,d,y_i);
  Real z_dT( unit = "1/K") = Molar.calc_z_dT_p(T,d,y_i);
  Real B( unit = "") = Molar.calc_B(T,d,y_i);
  Real B_dT( unit = "1/K") = Molar.calc_B_dT_p(T,d,y_i);
  SI.MolarVolume b = Molar.calc_b(y_i);
  Real a_alpha(unit="N.m4/mol2") = Molar.calc_a(T,y_i);
  Real a_alpha_dT(unit="N.m4/(mol2.K)") = Molar.calc_a_dT(T,y_i);
  Real a_alpha_d2T(unit="N.m4/(mol2.K2)") = Molar.calc_a_d2T(T,y_i);

algorithm
  cp := cp_ig + Modelica.Constants.R*(T*z_dT+z-1) + (T*a_alpha_dT-a_alpha)/(2*sqrt(2)*b)*((z_dT+(1+sqrt(2))*B_dT)/(z+(1+sqrt(2))*B) - (z_dT+(1-sqrt(2))*B_dT)/(z+(1-sqrt(2))*B)) + T*a_alpha_d2T/(2*sqrt(2)*b)*log((z+(1+sqrt(2))*B)/(z+(1-sqrt(2))*B));

end calc_cp;
