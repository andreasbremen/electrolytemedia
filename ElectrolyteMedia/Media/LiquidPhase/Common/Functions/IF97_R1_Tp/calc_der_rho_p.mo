within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.IF97_R1_Tp;
function calc_der_rho_p
  "Derivative w.r.t. p of density as function of pressure and temperaturefor region 1 of IF97"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  output Real rhop;
// protected
//   Real pi(unit="1") = Reduced.calc_pi(p) "Dimensionless pressure";
//   Real tau(unit="1") = Reduced.calc_tau(T) "Dimensionless temperature";
//   Real gpi(unit="1") = Reduced.calc_gpi(T, p) "Derivative of g w.r.t. pi";
//   Real gpipi(unit="1") = Reduced.calc_gpipi(T, p)
//     "2nd derivative of g w.r.t. pi";
//
// algorithm
//   rhop :=-1/(IF97.RH2O*T*gpi*gpi)*gpipi;
protected
  GibbsDerivs g = Reduced.GibbsDerivs(T,p);
  ThermoProperties_pT pro;
algorithm
  pro :=Modelica.Media.Common.ThermoFluidSpecial.gibbsToProps_pT(g);
  rhop :=pro.ddpT;
  annotation(smoothOrder = 5);
end calc_der_rho_p;
