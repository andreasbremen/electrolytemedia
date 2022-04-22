within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson.Molar;
function calc_lnphi
  "Calculates natural logarithm of fugacity coefficients of gas species with Peng Robinson EOS"

  input SI.Temperature T;
  input SI.Density d;
  input SI.MoleFraction[nGfun] y;

  output Real[nGfun] lnphi;
protected
  Real[nGfun] B_i(unit="") = Molar.calc_B_i(y);
  Real z(unit="") = Molar.calc_z(
    T,
    d,
    y);
  Real B(unit="") = Molar.calc_B(
    T,
    d,
    y);
  Real A(unit="") = Molar.calc_A(
    T,
    d,
    y);
  Real[nGfun] a_i(unit="N.m4/mol2") = Molar.calc_a_i_sum(T, y);
  Real a(unit="N.m4/mol2") = Molar.calc_a(T, y);
algorithm
  lnphi := B_i*(z - 1) - ones(nGfun)*log(
     z - B) + A/(2*sqrt(2)*B)*(B_i - 2*a_i/(a))*log((z + (1 + sqrt(2))*B)/(z -
    (sqrt(2) - 1)*B));
end calc_lnphi;
