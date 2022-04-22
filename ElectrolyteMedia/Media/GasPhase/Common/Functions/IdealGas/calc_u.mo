within ElectrolyteMedia.Media.GasPhase.Common.Functions.IdealGas;
function calc_u
  "Calculates specific internal energy of gas phase with ideal gas model"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nGfun] X;

  output SI.SpecificInternalEnergy u;
protected
  SI.SpecificEnthalpy[nGfun] hi=calc_h_i(T);
algorithm
  u :=X*hi - sum(datafun[i].R*X[i] for i in 1:nGfun)*T;

end calc_u;
