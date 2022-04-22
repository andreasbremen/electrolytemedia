within ElectrolyteMedia.Media.GasPhase.Common.Functions.IdealGas;
function calc_h
  "Calculates specific enthalpy of gas phase with ideal gas model"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nGfun] X;

  output SI.SpecificEnthalpy h;
protected
  SI.SpecificEnthalpy[nGfun] hi=calc_h_i(T);
algorithm
  h :=X*hi;

end calc_h;
