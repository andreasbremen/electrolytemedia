within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.SolventFunctions;
function calc_2der_g_TT
  "Second derivative of solvent function g w.r.t. T: d2gdT by Johnson1992, T and p dependent"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  output Real d2gdT;
protected
  Modelica.SIunits.Density rho = IF97_R1_Tp.calc_rho(T, p);
  Real rhoT=IF97_R1_Tp.calc_der_rho_T(T, p);
  Real rhoTT=IF97_R1_Tp.calc_2der_rho_TT(T, p);
  Modelica.SIunits.Temp_C T_C = T - 273.15;
  Real p_bar = p * 1e-5;
  Real rho_hat;
  Real a_g;
  Real b_g;
  Real da_g;
  Real d2a_g;
  Real db_g;
  Real d2b_g;
  Real g1;
  Real dg1dT;
  Real d2g1dT;
  Real d2ftdT;
  Real fp;
  Real d2fdT;
algorithm
//   if T_C > 1000 or p > 5000e5 or rho > 1000 or rho < 350 or T_C > 350 and T_C < 400 and p_bar < 500 then
//     Modelica.Utilities.Streams.error("Error: Model constraint violation. g function from Johnson et al. (1992) not valid for specified T and p. Provide different values or give better starting values.");
//   end if;
  if rho < 1000 then
    rho_hat := rho / Solvent.rho_0;
    a_g := Solvent.a[1] + Solvent.a[2] * T_C + Solvent.a[3] * T_C ^ 2;
    b_g := Solvent.b[1] + Solvent.b[2] * T_C + Solvent.b[3] * T_C ^ 2;
    da_g := Solvent.a[2] + 2 * Solvent.a[3] * T_C;
    db_g := Solvent.b[2] + 2 * Solvent.b[3] * T_C;
    d2a_g := 2 * Solvent.a[3];
    d2b_g := 2 * Solvent.b[3];
    g1 := (1 - rho_hat) ^ b_g;
    dg1dT :=g1*(db_g*log(1 - rho_hat) - b_g*rhoT/(Solvent.rho_0*(1 - rho_hat)));
    d2g1dT :=dg1dT*(db_g*log(1 - rho_hat) - b_g*rhoT/(Solvent.rho_0 - rho)) + g1*(d2b_g*
      log(1 - rho_hat) - 2*db_g*rhoT/(Solvent.rho_0 - rho) - b_g*rhoTT/(Solvent.rho_0 - rho) -
      b_g*rhoT^2/(Solvent.rho_0^2*(1 - rho_hat)^2));
    if T_C < 155 or T_C > 355 or p_bar > 1000 then
      d2ftdT := 0;
      fp := 0;
      d2fdT := 0;
    else
      fp := Solvent.a_f[2] * (1000 - p_bar) ^ 3 + Solvent.a_f[3] * (1000 - p_bar) ^ 4;
      d2ftdT := 3.8 * 4.8 / 300 ^ 2 * ((T_C - 155) / 300) ^ 2.8 + 15 * 16 / 300 ^ 2 * Solvent.a_f[1] * ((T_C - 155) / 300) ^ 14;
      d2fdT := d2ftdT * fp;
    end if;
    d2gdT := g1 * d2a_g + 2 * da_g * dg1dT + a_g * d2g1dT - d2fdT;
  else
    d2gdT := 0;
  end if;
  annotation(smoothOrder=20);
end calc_2der_g_TT;
