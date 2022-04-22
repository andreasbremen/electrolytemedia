within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel.Reduced.Additional.d2g_dtaudpi;
function calc_d2Aln_dtaudpi "Helper function to calculate Gibbs derivative"
  input SI.Temperature T;
  input SI.Pressure p;
  output Real dAln_dT;

protected
  SI.Density rho = IF97_R1_Tp.calc_rho(T,p);
  Real eps=BornFunctions.calc_eps(T,p);
  Real e=Modelica.Constants.F/Modelica.Constants.N_A "electronic charge 1.602176621e-19";
  Real tau=DebyeHueckel.Reduced.calc_tau(T);
  Real drho_dT=IF97_R1_Tp.calc_der_rho_T(T, p);
  Real d2rho_dTdp=IF97_R1_Tp.calc_2der_rho_pT(T, p);
  Real deps_dT=BornFunctions.calc_der_eps_T(T, p);
  Real d2eps_dTdp=BornFunctions.calc_2der_eps_pT(T, p);
  Real drho_dp=IF97_R1_Tp.calc_der_rho_p(T, p);
  Real deps_dp=BornFunctions.calc_der_eps_p(T, p);

algorithm

  dAln_dT := 1/3*(2*Modelica.Constants.pi*Modelica.Constants.N_A)^0.5*(
    e^2/(4*Modelica.Constants.pi*Modelica.Constants.epsilon_0*Modelica.Constants.k
    *Tref))^1.5*(0.5*((-0.5)*rho^(-1.5)*pref
    *drho_dp*(-Tref/tau^2)*drho_dT*(tau/eps)^1.5 + rho^(
    -0.5)*(-Tref/tau^2)*pref*d2rho_dTdp*
    (tau/eps)^1.5 + rho^(-0.5)*(-Tref/tau^2)*drho_dT*1.5
    *(tau/eps)^0.5*tau/eps^2*(-pref*deps_dp)) + 1.5*(0.5
    *rho^(-0.5)*pref*drho_dp*(tau/eps)^0.5*1/eps^2*(eps +
    Tref/tau*deps_dT) + rho^0.5*0.5*(tau/eps)^(-0.5)*
    tau/eps^2*(-pref*deps_dp)*1/eps^2*(eps + Tref
    /tau*deps_dT) + rho^0.5*(tau/eps)^0.5*(-2)*eps^(-3)*pref
    *deps_dp*(eps + Tref/tau*deps_dT) + rho^0.5*(tau/
    eps)^0.5*1/eps^2*(pref*deps_dp + Tref
    /tau*pref*d2eps_dTdp)));

end calc_d2Aln_dtaudpi;
