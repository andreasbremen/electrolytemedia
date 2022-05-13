within ElectrolyteMedia.Media.SolidLiquidPhase.Common.Functions;
function calc_Yfull "calculates total mole fraction vector"
  input SI.MassFraction [nLfunfun+nSfunfun] X;
  output SI.MoleFraction[nLfunfun+nSfunfun] Y;
protected
  Real[nLfunfun+nSfunfun] X_(unit="mol/kg");
algorithm
 for i in nSfunfun+1:nSfunfun+nLfunfun-1 loop
    X_[i] :=X[i]/dataLfunfun[i-nSfunfun].MM;
 end for;
 for i in 1:nSfunfun loop
    X_[i] :=X[i]/ dataSfunfun[i].MM;//0.05844;//
 end for;
  X_[nSfunfun+nLfunfun] :=X[nSfunfun+nLfunfun]/IF97.MH2O;
  Y :=X_/sum(X_);
end calc_Yfull;
