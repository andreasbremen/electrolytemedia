within ElectrolyteMedia.Media.GasLiquidPhase.Common.Functions;
function calc_X_M "calculates mass fraction vector of liquid species"
  input SI.MassFraction[:] Y;
  input SI.MolarMass[size(Y,1)] M;
  output SI.MoleFraction[size(Y,1)] X;
protected
  Real[size(Y,1)] Y_(unit="kg/mol");
algorithm
  for i in 1:size(Y,1) loop
    Y_[i] :=Y[i]*M[i];
  end for;
  X :=Y_/sum(Y_);

  annotation(smoothOrder=5);
end calc_X_M;
