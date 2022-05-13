within ElectrolyteMedia.Media.SolidLiquidPhase.Common.Functions;
function calc_Xl "calculates mass fraction vector of liquid phase"
  input SI.MassFraction[nSfunfun+nLfunfun] X;
  output SI.MassFraction[nLfunfun] Xl;
algorithm
  if sum(X[1+nSfunfun:nSfunfun+nLfunfun]) > 0 then
    Xl :=X[1+nSfunfun:nSfunfun+nLfunfun]/sum(X[1+nSfunfun:nSfunfun+nLfunfun]);
  else
    Xl :=fill(1/nLfunfun, nLfunfun);
  end if;

  annotation(smoothOrder=5);
end calc_Xl;
