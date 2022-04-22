within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.IF97_R1_Tp;
function calc_2der_rho_pT
  "Derivative w.r.t. p and T of density as function of pressure and temperaturefor region 1 of IF97"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  output Real rhopT;
protected
  Real pi(unit="1") = Reduced.calc_pi(p) "Dimensionless pressure";
  Real tau(unit="1") = Reduced.calc_tau(T) "Dimensionless temperature";
  Real gpi(unit="1") = Reduced.calc_der_g_pi(T, p) "Derivative of g w.r.t. pi";
  Real gpipi(unit="1") = Reduced.calc_2der_g_pipi(T, p)
    "2nd derivative of g w.r.t. pi";
  Real gtaupipi(unit="1") = Reduced.calc_3der_g_taupipi(T, p);
  Real gtaupi(unit="1") = Reduced.calc_2der_g_taupi(T, p);

algorithm

 rhopT :=(gpipi*(-2*tau*gtaupi+gpi)+tau*gpi*gtaupipi)/(IF97.RH2O*T^2*gpi*gpi*gpi);

  annotation(smoothOrder = 5);
end calc_2der_rho_pT;
