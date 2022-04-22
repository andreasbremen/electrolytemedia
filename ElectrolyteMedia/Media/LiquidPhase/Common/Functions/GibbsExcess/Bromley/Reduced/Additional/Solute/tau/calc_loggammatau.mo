within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Bromley.Reduced.Additional.Solute.tau;
function calc_loggammatau "Helper function to calculate Gibbs derivative"

  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nLfun] X;

  output Real[nLfun-1] log_gammatau;

protected
  Real Atau = DebyeHueckel.Reduced.Additional.dg_dtau.calc_dAlog_dtau(T,p);
  Real I=calc_I(X);
  Real[nLfun - 1] X_c=Bromley.Solute.calc_X_c(X);
  Real[nLfun - 1] X_a=Bromley.Solute.calc_X_a(X);
  Real[nLfun - 1] Y_ctau=calc_Y_ctau(T,p,X);
  Real[nLfun - 1] Y_atau=calc_Y_atau(T,p,X);
algorithm

  for i in 1:nLfun-1 loop
    if datafun[i].z > 0 then
      log_gammatau[i] :=Atau*I^0.5/(1 + I^0.5)*(-datafun[i].z^2 + X_c[i]) + Y_ctau[i];
    elseif datafun[i].z < 0 then
      log_gammatau[i] :=Atau*I^0.5/(1 + I^0.5)*(-datafun[i].z^2 + X_a[i]) + Y_atau[i];
    else
      log_gammatau[i] := 0.074 * I;
    end if;
  end for;

end calc_loggammatau;
