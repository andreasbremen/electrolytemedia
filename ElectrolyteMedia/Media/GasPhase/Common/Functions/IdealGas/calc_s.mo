within ElectrolyteMedia.Media.GasPhase.Common.Functions.IdealGas;
function calc_s
  "Calculates specific entropy of gas phase with ideal gas model"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nGfun] X;

  output SI.SpecificEntropy s;
protected
  SI.SpecificEntropy[nGfun] si=calc_s_i(T);
  SI.MoleFraction[nGfun] Y = calc_Y(X);
algorithm
  s:=X*si - sum(X[i]*Modelica.Constants.R/datafun[i].MM*
      (if X[i]<Modelica.Constants.eps then Y[i] else
      Modelica.Math.log(Y[i]*p/pref)) for i in 1:nGfun);

end calc_s;
