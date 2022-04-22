within ElectrolyteMedia.Media.SolidPhase.Common.Functions.Molar;
function calc_v_i "Calculates molar volume of mineral species"

  input SI.Temperature T;
  input SI.Pressure p;

  output SI.MolarVolume v[nSfun];
protected
  Real a[nSfun] = calc_a();
  Real b[nSfun] = calc_b();
  Real c[nSfun] = calc_c();
  SI.Pressure p_th[nSfun] = calc_p_th(T);
algorithm
  for i in 1:nSfun loop
    if datafun[:].k_0[i] == 0 then
      v[i] := datafun[:].V_ref[i];
    else
      v[i] := datafun[:].V_ref[i] * (1 - a[i] * (1 - (1 + b[i] * (p - p_th[i])) ^ (-c[i])));
    end if;
  end for;
  annotation(smoothOrder(normallyConstant = datafun)=5);
end calc_v_i;
