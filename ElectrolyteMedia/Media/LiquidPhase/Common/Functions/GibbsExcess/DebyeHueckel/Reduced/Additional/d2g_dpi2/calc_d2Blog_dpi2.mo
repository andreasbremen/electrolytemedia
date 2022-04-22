within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel.Reduced.Additional.d2g_dpi2;
function calc_d2Blog_dpi2 "Helper function to calculate Gibbs derivative"
  input SI.Temperature T;
  input SI.Pressure p;
  output Real d2Blog_dpi2;

protected
  SI.Density rho = IF97_R1_Tp.calc_rho(T,p);
  Real eps = BornFunctions.calc_eps(T,p);
  Real e=Modelica.Constants.F/Modelica.Constants.N_A
    "electronic charge 1.602176621e-19";
  Real drho_dp=IF97_R1_Tp.calc_der_rho_p(T, p);
  Real deps_dp=BornFunctions.calc_der_eps_p(T, p);
  Real d2rho_dp2=IF97_R1_Tp.calc_2der_rho_pp(T, p);
  Real d2eps_dp2=BornFunctions.calc_der_eps_pp(T, p);

algorithm
  d2Blog_dpi2 := (2*e^2*Modelica.Constants.N_A/(Modelica.Constants.epsilon_0
    *Modelica.Constants.k*T))^0.5*0.5*((-0.5)*(rho/eps)^(-1.5)*1/eps^4*(
    pref*drho_dp*eps - rho*pref*deps_dp)
    ^2 + (rho/eps)^(-0.5)*(-2)*eps^(-3)*pref*deps_dp*(
    pref*drho_dp*eps - rho*pref*deps_dp)
     + (rho/eps)^(-0.5)*1/eps^2*(pref^2*d2rho_dp2*eps -
    rho*pref^2*d2eps_dp2));

end calc_d2Blog_dpi2;
