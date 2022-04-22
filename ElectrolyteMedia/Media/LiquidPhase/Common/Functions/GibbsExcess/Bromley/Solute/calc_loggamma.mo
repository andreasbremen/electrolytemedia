within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Bromley.Solute;
function calc_loggamma
  "Calculates decadic logarithm of activity coefficient in molality base of aqueous species from extended Bromley equation, neutral species with b*I"

  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nLfun] X;

  output Real[nLfun-1] log_gamma;
protected
  Real DH=calc_DH_contribution(T,p,X);
  Real[nLfun-1] X_c = calc_X_c(X);
  Real[nLfun-1] X_a = calc_X_a(X);
  Real[nLfun-1] Y_c = calc_Y_c(T,p,X);
  Real[nLfun-1] Y_a = calc_Y_a(T,p,X);
  Real I=calc_I(X);
algorithm
  for i in 1:nLfun-1 loop
    if datafun[i].z > 0 then
      log_gamma[i] :=DH*(-datafun[i].z^2 + X_c[i]) + Y_c[i];
    elseif datafun[i].z < 0 then
      log_gamma[i] :=DH*(-datafun[i].z^2 + X_a[i]) + Y_a[i];
    else
      log_gamma[i] := 0.074 * I;
    end if;
  end for;
end calc_loggamma;
