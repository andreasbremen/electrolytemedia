within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Bromley.Reduced.Additional.Solute.taupi;
function calc_loggamma_ik_0taupi "Helper function to calculate Gibbs derivative"

  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nLfun] X;

  output Real[nLfun-1,nLfun-1] log_gammataupi;
protected
  Real Ataupi = DebyeHueckel.Reduced.Additional.d2g_dtaudpi.calc_d2Alog_dtaudpi(T,p);
  Real I=calc_I(X);
algorithm
  for i in 1:nLfun-1 loop
    for k in 1:nLfun-1 loop
      if datafun[i].z * datafun[k].z < 0 then
        log_gammataupi[i,k] := -Ataupi * abs(datafun[i].z* datafun[k].z) * I ^ 0.5 / (1 + I ^ 0.5);
      else
        log_gammataupi[i,k] :=0;
      end if;
    end for;
  end for;
end calc_loggamma_ik_0taupi;
