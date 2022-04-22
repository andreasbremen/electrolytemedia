within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer;
function calc_loggamma
  "Calculates decadic logarithm of activity coefficient in molality base of liquid species from extended Debye Hückel equation, neutral species with b*I"

  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nLfun] X;
  output Real[nLfun] loggamma;
protected
  SI.MoleFraction[nLfun] Y = calc_Y(X);
algorithm
  loggamma[1:nLfun-1] :=Solute.calc_loggamma(T,p,X);
  loggamma[nLfun] := Solvent.calc_loggamma(T,p,X);
  annotation(smoothOrder=20);
end calc_loggamma;
