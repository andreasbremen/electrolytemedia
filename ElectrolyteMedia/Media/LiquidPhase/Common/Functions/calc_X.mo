within ElectrolyteMedia.Media.LiquidPhase.Common.Functions;
function calc_X "Calculates mass fraction vector of liquid species from molar fraction vector"
  input SI.MassFraction[:] Y;
  output SI.MoleFraction[size(Y,1)] X;
protected
  Real[size(Y,1)] Y_(unit="kg/mol");
algorithm
  for i in 1:size(Y,1)-1 loop
    Y_[i] :=Y[i]*datafun[i].MM;
  end for;
  Y_[size(Y,1)] :=Y[size(Y,1)]*IF97.MH2O;
  X :=Y_/sum(Y_);

  annotation(smoothOrder=5);
end calc_X;
