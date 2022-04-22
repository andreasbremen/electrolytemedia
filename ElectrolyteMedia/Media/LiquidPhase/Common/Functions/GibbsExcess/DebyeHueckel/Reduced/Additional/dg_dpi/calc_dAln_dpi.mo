within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel.Reduced.Additional.dg_dpi;
function calc_dAln_dpi "Helper function to calculate Gibbs derivative"
  input SI.Temperature T;
  input SI.Pressure p;
  output Real dAln_dpi;

protected
  SI.Density rho = IF97_R1_Tp.calc_rho(T,p);
  Real eps=BornFunctions.calc_eps(T,p);
  Real e=Modelica.Constants.F/Modelica.Constants.N_A "electronic charge 1.602176621e-19";
  Real drho_dp=IF97_R1_Tp.calc_der_rho_p(T, p);
  Real deps_dp=BornFunctions.calc_der_eps_p(T, p);

algorithm

  dAln_dpi := 1/3*(2*Modelica.Constants.pi*Modelica.Constants.N_A)^0.5*(e^2/(4*
    Modelica.Constants.pi*Modelica.Constants.epsilon_0*Modelica.Constants.k*T))^
    1.5*(0.5*rho^(-0.5)*pref*drho_dp*eps^(-1.5) + rho^0.5*(-1.5)
    *eps^(-2.5)*pref*deps_dp);

end calc_dAln_dpi;
