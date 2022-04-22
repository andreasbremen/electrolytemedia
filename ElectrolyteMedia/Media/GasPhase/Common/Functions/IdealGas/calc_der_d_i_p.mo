within ElectrolyteMedia.Media.GasPhase.Common.Functions.IdealGas;
function calc_der_d_i_p
  "Calculates pressure derivative of density of ideal gas mixture at T and p"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nGfun] X;

  output SI.DerDensityByPressure ddpT;

algorithm
  ddpT :=1/(sum(X[i]*Modelica.Constants.R/datafun[i].MM for i in 1:nGfun)*T);
end calc_der_d_i_p;
