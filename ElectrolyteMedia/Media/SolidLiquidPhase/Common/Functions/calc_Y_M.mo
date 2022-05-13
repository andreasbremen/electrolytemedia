within ElectrolyteMedia.Media.SolidLiquidPhase.Common.Functions;
function calc_Y_M "calculates mole fraction vector of liquid species"
  input SI.MassFraction[:] X;
  input SI.MolarMass[size(X,1)] M;
  output SI.MoleFraction[size(X,1)] Y;
protected
  Real[size(X,1)] X_(unit="mol/kg");
algorithm
  for i in 1:size(X,1) loop
    X_[i] :=X[i]/M[i];
  end for;
  Y :=X_/sum(X_);

  annotation(smoothOrder=5);
end calc_Y_M;
