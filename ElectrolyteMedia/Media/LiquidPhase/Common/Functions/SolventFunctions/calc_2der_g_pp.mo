within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.SolventFunctions;
function calc_2der_g_pp
  "Derivative of solvent function g w.r.t. p: dgdp with IF97, T and p dependent"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  output Real dgd2p;
protected
  Modelica.SIunits.Density rho = IF97_R1_Tp.calc_rho(T, p);
  Real rhodp=IF97_R1_Tp.calc_der_rho_p(T, p);
  Real rhod2p=IF97_R1_Tp.calc_2der_rho_pp(T, p);
  Modelica.SIunits.Temp_C T_C = T - 273.15;
  Real p_bar = p * 1e-5;
  Real rho_hat;
  Real a_g;
  Real b_g;
  Real g1;
  Real f;
  Real ft;
  Real fp;
  Real dfpdp;
  Real dfdp;

  Real dg1d2p;
  Real dfpd2p;
  Real dfd2p;

  Real g;
algorithm
//   if T_C > 1000 or p > 5000e5 or rho > 1000 or rho < 350 or T_C > 350 and T_C < 400 and p_bar < 500 then
//     Modelica.Utilities.Streams.error("Error: Model constraint violation. g function from Johnson et al. (1992) not valid for specified T and p. Provide different values or give better starting values.");
//   end if;
  if rho < 1000 then
    rho_hat := rho / Solvent.rho_0;
    a_g := Solvent.a[1] + Solvent.a[2] * T_C + Solvent.a[3] * T_C ^ 2;
    b_g := Solvent.b[1] + Solvent.b[2] * T_C + Solvent.b[3] * T_C ^ 2;
    g1 := (1 - rho_hat) ^ b_g;
//     dg1d2p :=b_g*(1 - rho_hat)^(b_g - 2)*((b_g - 1)*(1/Solvent.rho_0*rhodp)^2 + (
//       rho_hat - 1)*(1 - Solvent.rho_0*rhod2p));
    dg1d2p :=(b_g*(1-rho/Solvent.rho_0)^b_g*((rho - Solvent.rho_0)*rhod2p + (b_g - 1)*rhodp^2))/(Solvent.rho_0 - rho)^2;
    if T_C < 155 or T_C > 355 or p_bar > 1000 then
      f := 0;
      ft := 0;
      fp := 0;
      dfpdp := 0;
      dfpd2p:=0;
      dfdp := 0;
      dfd2p :=0;
    else
      ft := ((T_C - 155) / 300) ^ 4.8 + Solvent.a_f[1] * ((T_C - 155) / 300) ^ 16;
      fp := Solvent.a_f[2] * (1000 - p_bar) ^ 3 + Solvent.a_f[3] * (1000 - p_bar) ^ 4;
      dfpd2p :=6*(1000 - p_bar)*(Solvent.a_f[2] + 2*Solvent.a_f[3]*(1000 - p_bar))*1e-10;
      f := ft * fp;
      dfpdp := ((-3 * Solvent.a_f[2] * (1000 - p_bar)^2) - 4 * Solvent.a_f[3] * (1000 - p_bar)^3) * 1e-5;
      dfdp := ft * dfpdp;
      dfd2p :=ft*dfpd2p;
    end if;
    g :=a_g*g1 - f;
    dgd2p := a_g * dg1d2p - dfd2p;
//     dgdp := (-g * rho_hat * beta * b_g / (1 - rho_hat)) - dfdp;
  else
    g := 0;
    dgd2p := 0;
  end if;
  annotation(smoothOrder=20);
end calc_2der_g_pp;
