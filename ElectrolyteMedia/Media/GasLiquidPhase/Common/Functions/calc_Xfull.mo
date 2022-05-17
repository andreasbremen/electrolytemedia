within ElectrolyteMedia.Media.GasLiquidPhase.Common.Functions;
function calc_Xfull
  "calculates full mass fraction vector from full mole fraction vector"
  input SI.MoleFraction[nLfunfun+nGfunfun] Y;
  output SI.MassFraction [nLfunfun+nGfunfun] X;

protected
  Real[nLfunfun+nGfunfun] Y_(unit="mol/kg");

algorithm
 for i in nGfunfun+1:nGfunfun+nLfunfun-1 loop
    Y_[i] :=Y[i]*dataLfunfun[i-nGfunfun].MM;
 end for;

 for i in 1:nGfunfun loop
    Y_[i] :=Y[i]* dataGfunfun[i].MM;//0.05844;//
 end for;
  Y_[nGfunfun+nLfunfun] :=Y[nGfunfun+nLfunfun]*IF97.MH2O;
  X :=Y_/sum(Y_);

end calc_Xfull;
