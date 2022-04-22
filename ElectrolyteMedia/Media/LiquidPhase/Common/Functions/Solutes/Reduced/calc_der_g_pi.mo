within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.Solutes.Reduced;
function calc_der_g_pi
  "Calculates derivative w.r.t. pi of Gibbs free energy at infinite dilution of aqueous species"

  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  output Real g[nLifun] "Molar Gibbs free energy at T, p";
protected
  Real tau=calc_tau(T);
  Real pi=calc_pi(p);
  Modelica.SIunits.MolarHeatCapacity R = Modelica.Constants.R;
  Real Z=BornFunctions.calc_Z(T, p);
  Real Zp=BornFunctions.calc_Z_dp(T, p);
  Real[nLifun] w=BornFunctions.calc_w(T,p);
  Real[nLifun] wp=BornFunctions.calc_der_w_p(T, p);
algorithm
  for i in 1:nLifun loop
    g[i] := datafun[i].a1*tau/R/Tref * pref + datafun[i].a2*tau/R/Tref * pref/(pi*pref+Psi) + tau^2 / (R*Tref*(Tref - Theta*tau)) * (datafun[i].a3 * pref + datafun[i].a4 * pref/(pi*pref+Psi)) - tau/R/Tref *pref*(wp[i]*(Z+1) + w[i]*Zp);
  end for;
//annotation(smoothOrder=5);
end calc_der_g_pi;
