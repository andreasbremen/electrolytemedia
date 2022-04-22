within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.IF97_R1_Tp;
function calc_2der_rho_TT
  "Calculate second derivative of density w.r.t. T from T and p with IF97"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;

  output Real rho_d2T;
protected
  Real tau=Reduced.calc_tau(T);
  Real pi=Reduced.calc_pi(p);
  Real gpi=Reduced.calc_der_g_pi(T, p);
  Real gtaupi=Reduced.calc_2der_g_taupi(T, p);
  Real gtautaupi=Reduced.calc_3der_g_tautaupi(T, p);
algorithm
  rho_d2T :=p*(2*tau^2*gtaupi^2 - 2*tau*gpi*gtaupi - 2*tau*gtaupi*gpi - tau^2*gtautaupi*gpi + 2*gpi^
    2)/(IF97.RH2O*T^3*gpi^3*pi);
  annotation(smoothOrder = 5);
end calc_2der_rho_TT;
