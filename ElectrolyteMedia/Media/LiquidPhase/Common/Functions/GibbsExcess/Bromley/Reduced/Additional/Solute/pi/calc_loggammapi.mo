within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Bromley.Reduced.Additional.Solute.pi;
function calc_loggammapi "Helper function to calculate Gibbs derivative"

  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nLfun] X;

  output Real[nLfun-1] log_gammapi;

protected
  Real Api = DebyeHueckel.Reduced.Additional.dg_dpi.calc_dAlog_dpi(T,p);
  Real I=calc_I(X);
  Real[nLfun - 1] X_c=Bromley.Solute.calc_X_c(X);
  Real[nLfun - 1] X_a=Bromley.Solute.calc_X_a(X);
  Real[nLfun-1] Y_cpi = calc_Y_cpi(T,p,X);
  Real[nLfun-1] Y_api = calc_Y_api(T,p,X);
algorithm

  for i in 1:nLfun-1 loop
    if datafun[i].z > 0 then
      log_gammapi[i] :=Api*I^0.5/(1 + I^0.5)*(-datafun[i].z^2 + X_c[i]) + Y_cpi[i];
    elseif datafun[i].z < 0 then
      log_gammapi[i] :=Api*I^0.5/(1 + I^0.5)*(-datafun[i].z^2 + X_a[i]) + Y_api[i];
    else
      log_gammapi[i] := 0.074 * I;
    end if;
  end for;

end calc_loggammapi;
