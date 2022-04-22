within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Bromley.Reduced.Additional.Solute.taupi;
function calc_loggammataupi "Helper function to calculate Gibbs derivative"

  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nLfun] X;

  output Real[nLfun-1] log_gammataupi;

protected
  Real Ataupi = DebyeHueckel.Reduced.Additional.d2g_dtaudpi.calc_d2Alog_dtaudpi(T,p);
  Real I=calc_I(X);
  Real[nLfun - 1] X_c=Bromley.Solute.calc_X_c(X);
  Real[nLfun - 1] X_a=Bromley.Solute.calc_X_a(X);
  Real[nLfun - 1] Y_ctaupi=calc_Y_ctaupi(T,p,X);
  Real[nLfun - 1] Y_ataupi=calc_Y_ataupi(T,p,X);
algorithm

  for i in 1:nLfun-1 loop
    if datafun[i].z > 0 then
      log_gammataupi[i] :=Ataupi*I^0.5/(1 + I^0.5)*(-datafun[i].z^2 + X_c[i]) + Y_ctaupi[i];
    elseif datafun[i].z < 0 then
      log_gammataupi[i] :=Ataupi*I^0.5/(1 + I^0.5)*(-datafun[i].z^2 + X_a[i]) + Y_ataupi[i];
    else
      log_gammataupi[i] := 0.074 * I;
    end if;
  end for;

end calc_loggammataupi;
