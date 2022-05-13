within ElectrolyteMedia.Media.SolidLiquidPhase.Common.Functions;
function calc_Xfull
  "calculates full mass fraction vector from full mole fraction vector"
  input SI.MoleFraction[nLfunfun+nSfunfun] Y;
  output SI.MassFraction [nLfunfun+nSfunfun] X;

protected
  Real[nLfunfun+nSfunfun] Y_(unit="mol/kg");

algorithm
 for i in nSfunfun+1:nSfunfun+nLfunfun-1 loop
    Y_[i] :=Y[i]*dataLfunfun[i-nSfunfun].MM;
 end for;

 for i in 1:nSfunfun loop
    Y_[i] :=Y[i]* dataSfunfun[i].MM;//0.05844;//
 end for;
  Y_[nSfunfun+nLfunfun] :=Y[nSfunfun+nLfunfun]*IF97.MH2O;
  X :=Y_/sum(Y_);

end calc_Xfull;
