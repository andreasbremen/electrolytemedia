within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Bromley;
function calc_gamma
  "Calculates activity coefficient in molality base of aqueous species from extended Debye Hückel equation, neutral species with b*I"

  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nLfun] X;
  output Real[nLfun] gamma;
protected
  Real[nLfun] loggamma;
algorithm

  loggamma:=calc_loggamma(T,p,X);
  for i in 1:nLfun loop
      gamma[i] :=10^loggamma[i];
  end for;
  annotation(smoothOrder=20);
end calc_gamma;
