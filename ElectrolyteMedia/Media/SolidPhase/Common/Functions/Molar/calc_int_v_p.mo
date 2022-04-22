within ElectrolyteMedia.Media.SolidPhase.Common.Functions.Molar;
function calc_int_v_p
  "Calculates integral of molar volume of mineral species over pressure"

  input SI.Temperature T;
  input SI.Pressure p;

  output Real vdp[nSfun];
protected
  parameter Real a[nSfun]=Molar.calc_a();
  parameter Real b[nSfun]=Molar.calc_b();
  parameter Real c[nSfun]=Molar.calc_c();
  SI.Pressure p_th[nSfun]=Molar.calc_p_th(T);
algorithm
  for i in 1:nSfun loop
    if datafun[:].k_0[i] == 0 then
      vdp[i] := datafun[:].V_ref[i] * (p - p_0);
    else
      vdp[i] := datafun[:].V_ref[i] * (p * (1 - a[i]) + a[i] * ((1 - b[i] * p_th[i]) ^ (1 - c[i]) - (1 + b[i] * (p - p_th[i])) ^ (1 - c[i])) / (b[i] * (c[i] - 1)));
    end if;
  end for;
  annotation(smoothOrder=20);
end calc_int_v_p;
