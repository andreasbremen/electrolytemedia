within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Bromley.Reduced.Additional.Solute.tautau;
function calc_loggammatautau "Helper function to calculate Gibbs derivative"

  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nLfun] X;

  output Real[nLfun-1] log_gammatautau;

protected
  Real Atautau = DebyeHueckel.Reduced.Additional.d2g_dtau2.calc_d2Alog_dtau2(T,p);
  Real I=calc_I(X);
  Real[nLfun - 1] X_c=Bromley.Solute.calc_X_c(X);
  Real[nLfun - 1] X_a=Bromley.Solute.calc_X_a(X);
  Real[nLfun - 1] Y_ctautau=calc_Y_ctautau(T,p,X);
  Real[nLfun - 1] Y_atautau=calc_Y_atautau(T,p,X);
algorithm

  for i in 1:nLfun-1 loop
    if datafun[i].z > 0 then
      log_gammatautau[i] :=Atautau*I^0.5/(1 + I^0.5)*(-datafun[i].z^2 + X_c[i]) + Y_ctautau[i];
    elseif datafun[i].z < 0 then
      log_gammatautau[i] :=Atautau*I^0.5/(1 + I^0.5)*(-datafun[i].z^2 + X_a[i]) + Y_atautau[i];
    else
      log_gammatautau[i] := 0.074 * I;
    end if;
  end for;

end calc_loggammatautau;
