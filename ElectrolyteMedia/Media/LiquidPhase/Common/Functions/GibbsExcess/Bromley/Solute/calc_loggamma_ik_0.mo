within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Bromley.Solute;
function calc_loggamma_ik_0
  "Calculates decadic logarithm of mean activity coefficient of single salt solution at ionic strength I"

  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nLfun] X;

  output Real[nLfun-1,nLfun-1] log_gamma;
protected
  Real A = DebyeHueckel.calc_A_log(T,p);
  Real I=calc_I(X);
algorithm
  for i in 1:nLfun-1 loop
    for k in 1:nLfun-1 loop
      if datafun[i].z * datafun[k].z < 0 then
        log_gamma[i,k] := -A * abs(datafun[i].z* datafun[k].z) * I ^ 0.5 / (1 + I ^ 0.5) + (0.06 + 0.6 * interactionfun.Bromley_ij[i,k]) * abs(datafun[i].z * datafun[k].z) * I / (1 + 1.5 ./ abs(datafun[i].z * datafun[k].z) * I) ^ 2 + interactionfun.Bromley_ij[i, k] * I;
      else
        log_gamma[i,k] :=0;
      end if;
    end for;
  end for;
end calc_loggamma_ik_0;
