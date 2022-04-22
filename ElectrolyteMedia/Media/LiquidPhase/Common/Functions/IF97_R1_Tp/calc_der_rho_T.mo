within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.IF97_R1_Tp;
function calc_der_rho_T
  "Calculate derivative of density w.r.t. T from T and p with IF97"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;

  output Real rho_dT;
// protected
//   Real tau=Reduced.calc_tau(T);
//   Real pi=Reduced.calc_pi(p);
//   Real gpi=Reduced.calc_gpi(T, p);
//   Real gtaupi=Reduced.calc_gtaupi(T, p);
// algorithm
//   rho_dT :=p/(IF97.RH2O*pi*T^2*gpi^2)*(tau*gtaupi - gpi);

protected
  GibbsDerivs g = Reduced.GibbsDerivs(T,p);
  ThermoProperties_pT pro;
algorithm
  pro :=Modelica.Media.Common.ThermoFluidSpecial.gibbsToProps_pT(g);
  rho_dT :=pro.ddTp;
  annotation(smoothOrder = 5);
end calc_der_rho_T;
