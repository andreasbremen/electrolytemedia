within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.SolventFunctions;
function calc_der_g_T
  "Derivative of solvent function g w.r.t. T: dgdT by Johnson1992, T and p dependent"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  output Real dgdT;
protected
  Modelica.SIunits.Density rho = IF97_R1_Tp.calc_rho(T, p);
  Real rhoT=IF97_R1_Tp.calc_der_rho_T(T, p);
  Modelica.SIunits.Temp_C T_C = T - 273.15;
  Real p_bar = p * 1e-5;
  Real rho_hat;
  Real a_g;
  Real b_g;
  Real da_g;
  Real db_g;
  Real g1;
  Real dg1dT;
  Real f;
  Real ft;
  Real dftdT;
  Real fp;
  Real dfdT;
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
    g1 := (1 - rho_hat) ^ b_g;
    dg1dT := g1 * (db_g * log(1 - rho_hat) - b_g * rhoT / (Solvent.rho_0*(1 - rho_hat)));
    if T_C < 155 or T_C > 355 or p_bar > 1000 then
      f := 0;
      ft := 0;
      dftdT := 0;
      fp := 0;
      dfdT := 0;
    else
      ft := ((T_C - 155) / 300) ^ 4.8 + Solvent.a_f[1] * ((T_C - 155) / 300) ^ 16;
      fp := Solvent.a_f[2] * (1000 - p_bar) ^ 3 + Solvent.a_f[3] * (1000 - p_bar) ^ 4;
      f := ft * fp;
      dftdT := 4.8 / 300 * ((T_C - 155) / 300) ^ 3.8 + 16 / 300 * Solvent.a_f[1] * ((T_C - 155) / 300) ^ 15;
      dfdT := dftdT * fp;
    end if;
    dgdT := a_g * dg1dT + da_g * g1 - dfdT;
  else
    dgdT := 0;
  end if;
  annotation(smoothOrder=20);
end calc_der_g_T;
