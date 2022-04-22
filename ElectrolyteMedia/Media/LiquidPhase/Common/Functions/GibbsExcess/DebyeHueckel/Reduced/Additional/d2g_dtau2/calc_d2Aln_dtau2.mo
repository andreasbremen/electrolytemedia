within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel.Reduced.Additional.d2g_dtau2;
function calc_d2Aln_dtau2 "Helper function to calculate Gibbs derivative"
  input SI.Temperature T;
  input SI.Pressure p;
  output Real d2Aln_dT2;

protected
  SI.Density rho = IF97_R1_Tp.calc_rho( T,p);
  Real eps=BornFunctions.calc_eps(T,p);
  Real e=Modelica.Constants.F/Modelica.Constants.N_A "electronic charge 1.602176621e-19";
  Real tau=DebyeHueckel.Reduced.calc_tau(T);
  Real drho_dT=IF97_R1_Tp.calc_der_rho_T(T, p);
  Real d2rho_dT2=IF97_R1_Tp.calc_2der_rho_TT(T, p);
  Real deps_dT=BornFunctions.calc_der_eps_T(T, p);
  Real d2eps_dT2=BornFunctions.calc_der_eps_TT(T, p);

algorithm

  d2Aln_dT2 := 1/3*(2*Modelica.Constants.pi*Modelica.Constants.N_A)^0.5*
    (e^2/(4*Modelica.Constants.pi*Modelica.Constants.epsilon_0*Modelica.Constants.k
    *Tref))^1.5*(0.5*((-0.5)*rho^(-1.5)*(-Tref
    /tau^2)^2*(drho_dT)^2*(tau/eps)^1.5 + rho^(-0.5)*2*Tref
    /tau^3*drho_dT*(tau/eps)^1.5 + rho^(-0.5)*(-Tref/
    tau^2)^2*d2rho_dT2*(tau/eps)^1.5 + rho^(-0.5)*(-Tref
    /tau^2)*drho_dT*1.5*(tau/eps)^0.5*1/eps^2*(eps + Tref
    /tau*deps_dT)) + 1.5*(0.5*rho^(-0.5)*(-Tref/tau^2)*
    drho_dT*(tau/eps)^0.5*((eps + Tref/tau*deps_dT)/eps^
    2) + rho^0.5*0.5*(tau/eps)^(-0.5)*(1/eps^2*(eps + Tref
    /tau*deps_dT))^2 + rho^0.5*(tau/eps)^0.5*1/eps^4*((-2*Tref
    /tau^2*deps_dT + Tref/tau*(-Tref/
    tau^2)*d2eps_dT2)*eps^2 - (eps + Tref/tau*deps_dT)*2
    *eps*(-Tref/tau^2)*deps_dT)));

end calc_d2Aln_dtau2;
