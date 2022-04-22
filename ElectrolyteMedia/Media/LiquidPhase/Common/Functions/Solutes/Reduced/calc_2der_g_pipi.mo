within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.Solutes.Reduced;
function calc_2der_g_pipi
  "Calculates second derivative w.r.t. pi of reduced Gibbs free energy at infinite dilution of aqueous species"

  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  output Real g[nLifun] "Molar Gibbs free energy at T, p";
protected
  Real tau=calc_tau(T);
  Real pi=calc_pi(p);
  Modelica.SIunits.MolarHeatCapacity R = Modelica.Constants.R;
  Real Z=BornFunctions.calc_Z(T, p);
  Real Zp=BornFunctions.calc_Z_dp(T, p);
  Real Zpp=BornFunctions.calc_Z_d2p(T, p);
  Real w[nLifun]=BornFunctions.calc_w(T,p);
  Real wp[nLifun]=BornFunctions.calc_der_w_p(T, p);
  Real wpp[nLifun]=BornFunctions.calc_2der_w_pp(T, p);
algorithm
  for i in 1:nLifun loop
    g[i] := -datafun[i].a2*tau/R/Tref * pref^2/(pi*pref+Psi)^2 - tau^2 / (R*Tref*(Tref - Theta*tau)) * datafun[i].a4 * pref^2/(pi*pref+Psi)^2 - tau/R/Tref *pref^2*(wpp[i]*(Z+1) + 2*wp[i]*Zp + w[i]*Zpp);
  end for;
//annotation(smoothOrder=5);
end calc_2der_g_pipi;
