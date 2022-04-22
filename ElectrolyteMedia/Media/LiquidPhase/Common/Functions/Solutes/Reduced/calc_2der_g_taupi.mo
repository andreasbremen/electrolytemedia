within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.Solutes.Reduced;
function calc_2der_g_taupi
  "Calculates second derivative w.r.t. tau and pi of Gibbs free energy at infinite dilution of aqueous species"

  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  output Real gtaupi[nLifun] "Molar Gibbs free energy at T, p";
protected
  Real tau=calc_tau(T);
  Real pi=calc_pi(p);
  Modelica.SIunits.MolarHeatCapacity R = Modelica.Constants.R;
  Real Z=BornFunctions.calc_Z(T, p);
  Real ZT=BornFunctions.calc_Z_dT(T, p);
  Real ZpT=BornFunctions.calc_Z_dpdT(T,p);
  Real Zp=BornFunctions.calc_Z_dp(T, p);
  Real w[nLifun]=BornFunctions.calc_w(T,p);
  Real wT[nLifun]=BornFunctions.calc_der_w_T(T, p);
  Real wpT[nLifun]=BornFunctions.calc_2der_w_pT(T, p);
  Real wp[nLifun]=BornFunctions.calc_der_w_p(T, p);
  Real Ttau = -Tref/tau^2;
  Real ppi = pref;
algorithm
  for i in 1:nLifun loop
    gtaupi[i] := datafun[i].a1/R/Tref * pref + datafun[i].a2/R/Tref * pref/(pi*pref+Psi) + tau*(2*Tref-Theta*tau)*pref / (R*Tref*(Tref - Theta*tau)^2) * (datafun[i].a3 + datafun[i].a4/(pi*pref+Psi)) - 1/R/Tref *ppi* (tau*Ttau*(Zp*wT[i]+(Z+1)*wpT[i]) + wp[i]*(tau*Ttau*ZT+Z+1) + w[i]*(tau*Ttau*ZpT+Zp));
  end for;
//annotation(smoothOrder=5);
end calc_2der_g_taupi;
