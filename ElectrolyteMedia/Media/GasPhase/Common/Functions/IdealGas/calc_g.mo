within ElectrolyteMedia.Media.GasPhase.Common.Functions.IdealGas;
function calc_g
  "Calculates specific Gibbs free energy of gas phase with ideal gas model"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nGfun] X;

  output SI.SpecificGibbsFreeEnergy g;
protected
  SI.SpecificGibbsFreeEnergy[nGfun] gi=calc_g_i(T);
  SI.MoleFraction[nGfun] Y = calc_Y(X);
algorithm
  g :=X*gi + T*sum(X[i]*Modelica.Constants.R/datafun[i].MM*
      (if X[i]<Modelica.Constants.eps then Y[i] else
      Modelica.Math.log(Y[i]*p/pref)) for i in 1:nGfun);

end calc_g;
