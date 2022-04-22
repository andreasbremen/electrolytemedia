within ElectrolyteMedia.Media.GasLiquidPhase.Common.Functions;
function calc_X "calculates mass fraction vector of liquid species"
  input SI.MoleFraction[nLfunfun] Y;
  output SI.MassFraction[nLfunfun] X;
protected
  Real[nGfunfun+nLfunfun] Y_(unit="kg/mol");
algorithm
  for i in 1:nGfunfun loop
    Y_[i] :=Y[i]*dataGfunfun[i].MM;
  end for;
  for i in 1:nLfunfun-1 loop
    Y_[i+nGfunfun] :=Y[i+nGfunfun]*dataLfunfun[i].MM;
  end for;
  Y_[nGfunfun+nLfunfun] := Y[nGfunfun+nLfunfun]*Common.IF97.MH2O;
  X :=Y_/sum(Y_);

  annotation(smoothOrder=5);
end calc_X;
