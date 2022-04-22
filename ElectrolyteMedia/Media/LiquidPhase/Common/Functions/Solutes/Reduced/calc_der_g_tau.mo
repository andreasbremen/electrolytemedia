within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.Solutes.Reduced;
function calc_der_g_tau
  "Calculates derivative w.r.t. tau of Gibbs free energy at infinite dilution of aqueous species"

  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  output Real gtau[nLifun] "Molar Gibbs free energy at T, p";
protected
  Real tau=calc_tau(T);
  Real pi=calc_pi(p);
  Modelica.SIunits.MolarHeatCapacity R = Modelica.Constants.R;
  Real Z=BornFunctions.calc_Z(T, p);
  Real ZT=BornFunctions.calc_Z_dT(T, p);
  Real w[nLifun]=BornFunctions.calc_w(T,p);
  Real wT[nLifun]=BornFunctions.calc_der_w_T(T, p);
  Real Ttau = -Tref/tau^2;
algorithm
  for i in 1:nLifun loop
      gtau[i] := datafun[i].H_ref/R/Tref + datafun[i].c1/R * (1/tau - 1) + datafun[i].a1/R/Tref * pref* (pi - 1) + datafun[i].a2/R/Tref * log((Psi + pi*pref) / (Psi + pref)) - datafun[i].c2/R/Tref * (Tref*(tau-1)/((Theta-Tref)*(Theta*tau-Tref)))   + tau*(2*Tref-Theta*tau) / (R*Tref*(Tref - Theta*tau)^2) * (datafun[i].a3 * pref*(pi - 1) + datafun[i].a4 * log((Psi + pi*pref) / (Psi + pref))) - 1/R/Tref * (tau*(Z+1)*Ttau*wT[i] + w[i]*(tau*Ttau*ZT+Z+1)) + datafun[i].w_ref/R/Tref * (Zref + 1) - datafun[i].w_ref/R * Yref;
  end for;
//annotation(smoothOrder=5);
end calc_der_g_tau;
