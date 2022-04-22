within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.BornFunctions;
function calc_der_w_p
  "Calculates derivative of omega function w.r.t. p: dwdp by Shock et al. 1992"

  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  output Real[nLifun] dwdp;
protected
  Real g = SolventFunctions.calc_g(T,p);
  Real dgdp=SolventFunctions.calc_der_g_p(T, p);
  Real[nLifun] z_abs;
  Real[nLifun] w_ref_cal;
  Real[nLifun] r_ref;
  Real[nLifun] r;
algorithm
  for i in 1:nLifun loop
    z_abs[i] :=abs(datafun[i].z);
    w_ref_cal[i] :=datafun[i].w_ref/4.184;
    if datafun[i].z == 0 or datafun[i].MM < 0.001008 then
      r_ref[i] := 0;
      r[i] := 0;
      dwdp[i] := 0;
    else
      r_ref[i] := datafun[i].z ^ 2 / (w_ref_cal[i] / eta + datafun[i].z / 3.082);
      r[i] := r_ref[i] + z_abs[i] * g;
      dwdp[i] := (-eta * (z_abs[i] ^ 3 / r[i] ^ 2 - datafun[i].z / (3.082 + g) ^ 2) * dgdp) * 4.184;
    end if;
  end for;
  annotation(smoothOrder=20);
end calc_der_w_p;
