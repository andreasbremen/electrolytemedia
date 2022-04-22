within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson.Molar;
function calc_A_dT_p
  "Calculates derivative w.r.t. T at const p of product of van der Waals attraction and acentric factor over squared RT in Peng Robinson EOS"

  input SI.Temperature T;
  input SI.Density d;
  input SI.MoleFraction[nGfun] y_i;

  output Real A_dT(unit="1/K");
protected
  SI.Pressure p=Molar.calc_p(T,d,y_i);
  Real a_alpha(unit="N.m4/mol2") = Molar.calc_a(T, y_i);
  Real a_alpha_dT(unit="N.m4/(mol2.K)") = Molar.calc_a_dT(T, y_i);
algorithm
  A_dT:=p/(Modelica.Constants.R*T)^2*(a_alpha_dT - 2*a_alpha/T);
end calc_A_dT_p;
