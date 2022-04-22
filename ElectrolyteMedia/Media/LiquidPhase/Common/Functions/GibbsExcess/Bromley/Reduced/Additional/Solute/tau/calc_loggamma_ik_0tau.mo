within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Bromley.Reduced.Additional.Solute.tau;
function calc_loggamma_ik_0tau "Helper function to calculate Gibbs derivative"

  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nLfun] X;

  output Real[nLfun-1,nLfun-1] log_gammatau;
protected
  Real Atau = DebyeHueckel.Reduced.Additional.dg_dtau.calc_dAlog_dtau(T,p);
  Real I=calc_I(X);
algorithm
  for i in 1:nLfun-1 loop
    for k in 1:nLfun-1 loop
      if datafun[i].z * datafun[k].z < 0 then
        log_gammatau[i,k] := -Atau * abs(datafun[i].z* datafun[k].z) * I ^ 0.5 / (1 + I ^ 0.5);
      else
        log_gammatau[i,k] :=0;
      end if;
    end for;
  end for;
end calc_loggamma_ik_0tau;
