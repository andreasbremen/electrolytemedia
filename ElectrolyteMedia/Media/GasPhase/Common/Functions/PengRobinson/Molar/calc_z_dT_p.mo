within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson.Molar;
function calc_z_dT_p
  "Calculates derivative w.r.t. T at const p of compressibility factor of gas mixture with Peng Robinson EOS"

  input SI.Temperature T;
  input SI.Density d;
  input SI.MoleFraction[nGfun] y_i;

  output Real z_dT(unit="");
protected
  SI.MolarMass MM = Molar.calc_MM(y_i);
  SI.MolarVolume v = MM/d;
  Real z=Molar.calc_z(
      T,
      d,
      y_i);
  Real A=Molar.calc_A(
      T,
      d,
      y_i);
  Real A_dT=Molar.calc_A_dT_p(
      T,
      d,
      y_i);
  Real B=Molar.calc_B(
      T,
      d,
      y_i);
  Real B_dT=Molar.calc_B_dT_p(
      T,
      d,
      y_i);
algorithm
  z_dT :=(A_dT*(B - z) + B_dT*(6*B*z + 2*z - 3*B^2 - 2*B + A - z^2))/(3*z^2 + 2*
    (B - 1)*z + (A - 2*B - 3*B^2));
end calc_z_dT_p;
