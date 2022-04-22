within ElectrolyteMedia.Media.GasLiquidPhase.Common.Functions;
function calc_Y "calculates mole fraction vector of liquid species"
  input SI.MassFraction[nGfunfun+nLfunfun] X;
  output SI.MoleFraction[nGfunfun+nLfunfun] Y;
protected
  Real[nGfunfun+nLfunfun] X_(unit="mol/kg");
algorithm
  for i in 1:nGfunfun loop
    X_[i] :=X[i]/dataGfunfun[i].MM;
  end for;
  for i in 1:nLfunfun-1 loop
    X_[i+nGfunfun] :=X[i+nGfunfun]/dataLfunfun[i].MM;
  end for;
  X_[nGfunfun+nLfunfun] := X[nGfunfun+nLfunfun]/Common.IF97.MH2O;
  Y :=X_/sum(X_);

  annotation(smoothOrder=5);
end calc_Y;
