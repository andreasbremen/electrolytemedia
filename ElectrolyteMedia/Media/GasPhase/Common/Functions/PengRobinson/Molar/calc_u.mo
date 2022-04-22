within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson.Molar;
function calc_u
  "Calculates internal energy of gas phase with Peng Robinson"
  input SI.Temperature T;
  input SI.Density d;
  input SI.MoleFraction[nGfun] y_i;

  output SI.MolarEnergy u;

protected
    SI.MolarMass MM = calc_MM(y_i);
    SI.MolarVolume v = MM/d;
  SI.MolarEnthalpy[nGfun] h_ig=IdealGas.Molar.calc_h_i(T);
    SI.MolarEnergy u_ig= y_i*h_ig - Modelica.Constants.R*T;
    Real a(unit="N.m4/mol2") = calc_a(T,y_i);
    Real a_dT(unit="N.m4/(mol2.K)")=calc_a_dT(T,y_i);
    SI.MolarVolume b=calc_b(y_i);
algorithm
  u := u_ig + (a-T*a_dT)/(2*sqrt(2)*b)*log((v+(1-sqrt(2))*b)/(v+(1+sqrt(2))*b));

end calc_u;
