within ElectrolyteMedia.Media.GasPhase.Common.Functions;
function calc_X
  "Calculates mass fraction vector from mole fraction vector"

  input SI.MassFraction[nGfun] X;
  output SI.MoleFraction[nGfun] Y;

protected
  Real[nGfun] Y_ = X.*datafun[:].MM;

algorithm
  Y :=Y_/sum(Y_);

end calc_X;
