within ElectrolyteMedia.Media.GasPhase.Common.Functions.IdealGas;
function calc_der_d_i_T
  "Calculates temperature derivative of density of ideal gas mixture at T and p"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nGfun] X;

  output SI.DerDensityByTemperature ddTp;

algorithm
  ddTp :=-p/(T*T*sum(X[i]*Modelica.Constants.R/datafun[i].MM for i in 1:nGfun));
end calc_der_d_i_T;
