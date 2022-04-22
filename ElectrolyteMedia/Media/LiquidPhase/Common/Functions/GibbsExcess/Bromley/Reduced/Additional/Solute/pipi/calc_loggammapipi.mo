within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Bromley.Reduced.Additional.Solute.pipi;
function calc_loggammapipi "Helper function to calculate Gibbs derivative"

  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nLfun] X;

  output Real[nLfun-1] log_gammapipi;

protected
  Real Apipi = DebyeHueckel.Reduced.Additional.d2g_dpi2.calc_d2Alog_dpi2(T,p);
  Real I=calc_I(X);
  Real[nLfun - 1] X_c=Bromley.Solute.calc_X_c(X);
  Real[nLfun - 1] X_a=Bromley.Solute.calc_X_a(X);
  Real[nLfun - 1] Y_cpipi=calc_Y_cpipi(T,p,X);
  Real[nLfun - 1] Y_apipi=calc_Y_apipi(T,p,X);
algorithm

  for i in 1:nLfun-1 loop
    if datafun[i].z > 0 then
      log_gammapipi[i] :=Apipi*I^0.5/(1 + I^0.5)*(-datafun[i].z^2 + X_c[i]) + Y_cpipi[i];
    elseif datafun[i].z < 0 then
      log_gammapipi[i] :=Apipi*I^0.5/(1 + I^0.5)*(-datafun[i].z^2 + X_a[i]) + Y_apipi[i];
    else
      log_gammapipi[i] := 0.074 * I;
    end if;
  end for;

end calc_loggammapipi;
