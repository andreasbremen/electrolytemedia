within ElectrolyteMedia.Media.SolidLiquidPhase.Common.Functions;
function calc_Xs "calculates mass fraction vector of solid phase"
  input SI.MassFraction[nSfunfun+nLfunfun] X;
  output SI.MassFraction[nSfunfun] Xs;
algorithm
  if sum(X[1:nSfunfun]) > 0 then
    Xs :=X[1:nSfunfun]/sum(X[1:nSfunfun]);
  else
    Xs :=fill(1/nSfunfun, nSfunfun);
  end if;

  annotation(smoothOrder=5);
end calc_Xs;
