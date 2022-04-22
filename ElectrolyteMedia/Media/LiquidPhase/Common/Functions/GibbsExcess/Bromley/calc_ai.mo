within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Bromley;
function calc_ai "Calculates activity in molality base with Bromley model"

  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nLfun] X;

  output Real[nLfun] ai;

protected
  Real[nLfun-1] loggamma = Solute.calc_loggamma(T,p,X);
  Real[nLfun - 1] molality=calc_mfromX(X);
  Real loggamma_s = Solvent.calc_loggamma(T,p,X);

algorithm

  for i in 1:nLfun-1 loop
    ai[i] := 10^loggamma[i] * molality[i];
  end for;

  ai[nLfun] := 10^loggamma_s/IF97.MH2O;

end calc_ai;
