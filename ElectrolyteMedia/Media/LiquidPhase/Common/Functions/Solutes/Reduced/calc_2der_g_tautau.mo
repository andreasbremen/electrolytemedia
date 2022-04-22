within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.Solutes.Reduced;
function calc_2der_g_tautau
  "Calculates second derivative w.r.t. tau of Gibbs free energy at infinite dilution of aqueous species"

  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  output Real gtautau[nLifun] "Molar Gibbs free energy at T, p";
protected
  Real tau=calc_tau(T);
  Real pi=calc_pi(p);
  Modelica.SIunits.MolarHeatCapacity R = Modelica.Constants.R;
  Real Z=BornFunctions.calc_Z(T, p);
  Real ZT=BornFunctions.calc_Z_dT(T, p);
  Real ZTT=BornFunctions.calc_Z_d2T(T, p);
  Real w[nLifun]=BornFunctions.calc_w(T,p);
  Real wT[nLifun]=BornFunctions.calc_der_w_T(T, p);
  Real wTT[nLifun]=BornFunctions.calc_2der_w_TT(T, p);
  Real Ttau = -Tref/tau^2;
algorithm
  for i in 1:nLifun loop
    gtautau[i] := -datafun[i].c1/R /tau^2 - datafun[i].c2/R/((Theta*tau-Tref)^2)   + 2*Tref / (R*(Tref - Theta*tau)^3) * (datafun[i].a3 * pref*(pi - 1) + datafun[i].a4 * log((Psi + pi*pref) / (Psi + pref))) - 1/R *(1/tau^2*wT[i]*(Z+1) - (Z+1)*wTT[i]*Ttau/tau - wT[i]*ZT*Ttau/tau + wT[i]*Ttau*(-1/tau*ZT+Z/Tref+1/Tref) + w[i]*(ZT/tau^2-1/tau*Ttau*ZTT+ZT*Ttau/Tref));
  end for;
  // gtautau := -datafun[i].c1/R /tau^2 - datafun[i].c2/R/((Theta*tau-Tref)^2)   + 2*Tref / (R*(Tref - Theta*tau)^3) * (datafun[i].a3 * pref*(pi - 1) + datafun[i].a4 * log((Psi + pi*pref) / (Psi + pref))) - 1/R *(1/tau^2*wT*(Z+1) - (Z+1)*wTT*Ttau - wT*ZT*Ttau + wT*Ttau*(-1/tau*ZT+Z+1) + w*(ZT/tau^2-1/tau*Ttau*ZTT+ZT*Ttau));
//annotation(smoothOrder(normallyConstant = datafun[i])=20);
end calc_2der_g_tautau;
