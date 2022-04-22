within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.IF97_R1_Tp;
function calc_2der_rho_pp
  "Second derivative w.r.t. p of density as function of pressure and temperature for region 1 of IF97"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  output Real rhopp;
protected
  Real pi(unit="1") = Reduced.calc_pi(p) "Dimensionless pressure";
  Real tau(unit="1") = Reduced.calc_tau(T) "Dimensionless temperature";
  Real gpi(unit="1") = Reduced.calc_der_g_pi(T, p) "Derivative of g w.r.t. pi";
  Real gpipi(unit="1") = Reduced.calc_2der_g_pipi(T, p)
    "2nd derivative of g w.r.t. pi";
  Real gpipipi(unit="1") = Reduced.calc_3der_g_pipipi(T, p)
    "3nd derivative of g w.r.t. pi";

algorithm

  rhopp :=pi/(IF97.RH2O*T*p*gpi^3)*(2*gpipi^2 - gpi*gpipipi);

  annotation(smoothOrder = 5);
end calc_2der_rho_pp;
