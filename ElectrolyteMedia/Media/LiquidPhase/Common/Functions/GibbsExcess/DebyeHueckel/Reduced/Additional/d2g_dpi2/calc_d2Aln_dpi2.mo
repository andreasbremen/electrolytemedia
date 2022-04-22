within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel.Reduced.Additional.d2g_dpi2;
function calc_d2Aln_dpi2 "Helper function to calculate Gibbs derivative"
  input SI.Temperature T;
  input SI.Pressure p;
  output Real d2Aln_dpi2;

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

  d2Aln_dpi2 := 1/3*(2*Modelica.Constants.pi*Modelica.Constants.N_A)^0.5
    *(e^2/(4*Modelica.Constants.pi*Modelica.Constants.epsilon_0*
    Modelica.Constants.k*T))^1.5*(0.5*(-0.5*rho^(-1.5)*pref
    ^2*(drho_dp)^2*eps^(-1.5) + rho^(-0.5)*pref^2*
    d2rho_dp2*eps^(-1.5) + rho^(-0.5)*pref^2*drho_dp*(-1.5)
    *eps^(-2.5)*deps_dp) - 1.5*(0.5*rho^(-0.5)*pref^2*
    drho_dp*eps^(-2.5)*deps_dp + rho^(0.5)*(-2.5)*eps^(-3.5)*pref
    ^2*deps_dp^2 + rho^(0.5)*eps^(-2.5)*pref^2*
    d2eps_dp2));

end calc_d2Aln_dpi2;
