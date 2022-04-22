within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.Solutes.Reduced;
function calc_g
  "Calculates reduced Gibbs free energy at infinite dilution of aqueous species"

  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  output Real[nLifun] g;
protected
  Real tau=calc_tau(T);
  Real pi=calc_pi(p);
  Modelica.SIunits.MolarHeatCapacity R = Modelica.Constants.R;
  Modelica.SIunits.Density rho = IF97_R1_Tp.calc_rho(T,p);
  Real Z=BornFunctions.calc_Z(T, p);
  Real[nLifun] w=BornFunctions.calc_w(T,p);
algorithm
  for i in 1:nLifun loop
    g[i] := datafun[i].G_ref*tau/R/Tref - datafun[i].S_ref/R*(1-tau) - datafun[i].c1/R * (log(1/tau) - 1 + tau) + datafun[i].a1*tau/R/Tref * pref* (pi - 1) + datafun[i].a2*tau/R/Tref * log((Psi + pi*pref) / (Psi + pref)) - datafun[i].c2/R/Tref * ((tau / (Tref - Theta*tau) - 1 / (Tref - Theta)) * (Theta*tau - Tref) / Theta - Tref / Theta ^ 2 * log((Tref - Theta*tau) / (Tref - Theta))) + tau^2 / (R*Tref*(Tref - Theta*tau)) * (datafun[i].a3 * pref*(pi - 1) + datafun[i].a4 * log((Psi + pi*pref) / (Psi + pref))) - w[i]*tau/R/Tref * (Z + 1) + datafun[i].w_ref*tau/R/Tref * (Zref + 1) + datafun[i].w_ref/R * Yref * (1 - tau);
  end for;
//annotation(smoothOrder=5);
end calc_g;
