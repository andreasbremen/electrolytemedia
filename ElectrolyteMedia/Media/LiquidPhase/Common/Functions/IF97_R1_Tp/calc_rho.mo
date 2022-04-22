within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.IF97_R1_Tp;
function calc_rho "Calculates density from temperature and pressure"
  input SI.Temperature T;
  input SI.Pressure p;
  output SI.Density rho;
// protected
//   Real pi=Reduced.calc_pi(p);
//   Real gpi=Reduced.calc_gpi(T, p);
// algorithm
//   rho :=p/(IF97.RH2O*T*pi*gpi);
protected
  GibbsDerivs g = Reduced.GibbsDerivs(T,p);
  ThermoProperties_pT pro;
algorithm
  pro :=Modelica.Media.Common.ThermoFluidSpecial.gibbsToProps_pT(g);
  rho :=pro.d;
  annotation(smoothOrder = 5);
end calc_rho;
