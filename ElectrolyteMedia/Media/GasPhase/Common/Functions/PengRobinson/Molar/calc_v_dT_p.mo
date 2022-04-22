within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson.Molar;
function calc_v_dT_p "Calculates derivative of molar volume w.r.t. T at constant p"
    input SI.Temperature T;
  input SI.Density d;
  input SI.MoleFraction[nGfun] y_i;

  output Real v_dT(unit="m3/(mol.K)");
protected
  SI.Pressure p = Molar.calc_p(T,d,y_i);
  Real z=Molar.calc_z(
      T,
      d,
      y_i);
  Real z_dT=Molar.calc_z_dT_p(
      T,
      d,
      y_i);
algorithm
  v_dT :=Modelica.Constants.R/p*(T*z_dT + z);

end calc_v_dT_p;
