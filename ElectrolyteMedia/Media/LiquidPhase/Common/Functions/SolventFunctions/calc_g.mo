within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.SolventFunctions;
function calc_g
  "Solvent function g and its derivatives dgdT, d2gdT, dgdp determined by Johnson1992, T and p dependent"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  output Real g;
protected
  Modelica.SIunits.Density rho = IF97_R1_Tp.calc_rho(T, p);
  Modelica.SIunits.Temp_C T_C = T - 273.15;
  Real p_bar = p * 1e-5;
  Real rho_hat;
  Real a_g;
  Real b_g;
  Real g1;
  Real f;
  Real ft;
  Real fp;
algorithm
//   if T_C > 1000 or p > 5000e5 or rho > 1000 or rho < 350 or (T_C > 350 and T_C < 400 and p_bar < 500) then
//     Modelica.Utilities.Streams.error("Error: Model constraint violation. g function from Johnson et al. (1992) not valid for specified T and p. Provide different values or give better starting values.");
//   end if;
  if rho < 1000 then
    rho_hat := rho / Solvent.rho_0;
    a_g := Solvent.a[1] + Solvent.a[2] * T_C + Solvent.a[3] * T_C ^ 2;
    b_g := Solvent.b[1] + Solvent.b[2] * T_C + Solvent.b[3] * T_C ^ 2;
    g1 := (1 - rho_hat) ^ b_g;
    if T_C < 155 or T_C > 355 or p_bar > 1000 then
      ft := 0;
      fp := 0;
      f := 0;
    else
      ft := ((T_C - 155) / 300) ^ 4.8 + Solvent.a_f[1] * ((T_C - 155) / 300) ^ 16;
      fp := Solvent.a_f[2] * (1000 - p_bar) ^ 3 + Solvent.a_f[3] * (1000 - p_bar) ^ 4;
      f := ft * fp;
    end if;
    g := a_g * g1 - f;
  else
    g := 0;
  end if;
  annotation(smoothOrder=20);
end calc_g;
