within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel.Reduced.Additional.dg_dtau;
function calc_dAln_dtau "Helper function to calculate Gibbs derivative"
  input SI.Temperature T;
  input SI.Pressure p;
  output Real dAln_dT;

protected
  SI.Density rho = IF97_R1_Tp.calc_rho(T,p);
  Real eps=BornFunctions.calc_eps(T,p);
  Real e=Modelica.Constants.F/Modelica.Constants.N_A "electronic charge 1.602176621e-19";
  Real tau=calc_tau(T);
  Real drho_dT=IF97_R1_Tp.calc_der_rho_T(T, p);
  Real deps_dT=BornFunctions.calc_der_eps_T(T, p);

algorithm

  dAln_dT := 1/3*(2*Modelica.Constants.pi*Modelica.Constants.N_A)^0.5*(
    e^2/(4*Modelica.Constants.pi*Modelica.Constants.epsilon_0*Modelica.Constants.k
    *Tref))^1.5*(0.5*rho^(-0.5)*(-Tref/
    tau^2)*drho_dT*(tau/eps)^1.5 + rho^0.5*1.5*(tau/eps)^0.5*((eps - (-
    Tref/tau)*deps_dT)/eps^2));

end calc_dAln_dtau;
