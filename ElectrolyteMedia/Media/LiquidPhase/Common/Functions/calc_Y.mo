within ElectrolyteMedia.Media.LiquidPhase.Common.Functions;
function calc_Y "Calculates mole fraction vector of liquid species"
  input SI.MassFraction[nLfun] X;
  output SI.MoleFraction[nLfun] Y;
protected
  Real[nLfun] X_(unit="mol/kg");
algorithm
  for i in 1:nLfun-1 loop
    X_[i] :=X[i]/datafun[i].MM;
  end for;
  X_[nLfun] :=X[nLfun]/IF97.MH2O;
  Y :=X_/sum(X_);

  annotation(smoothOrder=5);
end calc_Y;
