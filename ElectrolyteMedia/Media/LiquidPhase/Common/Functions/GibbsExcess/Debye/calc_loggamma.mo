within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Debye;
function calc_loggamma
  "Calculates decadic logarithm of activity coefficient in molality base of liquid species from Debye limiting law"

  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nLfun] X;
  output Real[nLfun] loggamma;
algorithm
  loggamma[1:nLfun-1] :=calc_loggamma_i(T,p,X);
  loggamma[nLfun] := Solvent.calc_loggamma(T,p,X);//log10(Y[nLfun]*MixtureLiquid.MH2O);
  annotation(smoothOrder=20);
end calc_loggamma;
