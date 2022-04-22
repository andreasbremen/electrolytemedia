within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel.Reduced.Additional.dg_dpi;
function calc_dBlog_dpi "Helper function to calculate Gibbs derivative"
  input SI.Temperature T;
  input SI.Pressure p;
  output Real dBlog_dpi;

protected
  SI.Density rho= IF97_R1_Tp.calc_rho(T,p);
  Real eps = BornFunctions.calc_eps(T,p);
  Real e=Modelica.Constants.F/Modelica.Constants.N_A
    "electronic charge 1.602176621e-19";
  Real drho_dp=IF97_R1_Tp.calc_der_rho_p(T, p);
  Real deps_dp=BornFunctions.calc_der_eps_p(T, p);
algorithm
  dBlog_dpi := (2*e^2*Modelica.Constants.N_A/(Modelica.Constants.epsilon_0*
    Modelica.Constants.k*T))^0.5*(0.5*(rho/eps)^(-0.5)*1/eps^2*(pref
    *drho_dp*eps - rho*pref*deps_dp));

end calc_dBlog_dpi;
