within ElectrolyteMedia.Media.LiquidPhase.Common.Functions;
function calc_Xfromm
  "Calculates mass fraction vector from molality vector of liquid species"
  input Real[nLfun-1] m;
  output SI.MassFraction[nLfun] X;
protected
  SI.MoleFraction[nLfun] Y = calc_Yfromm(m);
algorithm
  for i in 1:nLfun-1 loop
    X[i] :=Y[i]*datafun[i].MM/(Y[1:nLfun-1]*datafun[:].MM + Y[nLfun]*IF97.MH2O);
  end for;
  X[nLfun] :=1-sum(X[1:nLfun - 1]);

  annotation(smoothOrder=5);
end calc_Xfromm;
