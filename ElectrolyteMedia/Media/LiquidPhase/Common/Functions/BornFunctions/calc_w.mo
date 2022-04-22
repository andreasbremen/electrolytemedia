within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.BornFunctions;
function calc_w
  "Calculates omega function w determined by Shock et al. 1992"

  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  output Real[nLifun] w;
protected
  Real g = SolventFunctions.calc_g(T,p);
  Real[nLifun] z_abs;
  Real[nLifun] w_ref_cal;
  Real[nLifun] r_ref;
  Real[nLifun] r;
algorithm
  for i in 1:nLifun loop
    z_abs[i] :=datafun[i].z^2;
    w_ref_cal[i] :=datafun[i].w_ref/4.184;
    if datafun[i].z == 0 or datafun[i].MM < 0.001008 then
      r_ref[i] := 0;
      r[i] := 0;
      w[i] := w_ref_cal[i]*4.184;
    else
      r_ref[i] := datafun[i].z ^ 2 / (w_ref_cal[i] / eta + datafun[i].z / 3.082);
      r[i] := r_ref[i] + z_abs[i] * g;
      w[i] := eta * (datafun[i].z ^ 2 / r[i] - datafun[i].z / (3.082 + g)) *4.184;
    end if;
  end for;
  annotation(smoothOrder = 5);
end calc_w;
