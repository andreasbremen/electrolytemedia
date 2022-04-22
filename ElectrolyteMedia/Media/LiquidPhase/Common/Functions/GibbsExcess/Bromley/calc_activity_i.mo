within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Bromley;
function calc_activity_i "Calculates activity in molality base with Debye Hückel model"

input SI.Temperature T;
input SI.Pressure p;
input SI.MassFraction[nLfun] X;
output Real[nLfun-1] ai;

protected
Real[nLfun - 1] loggamma=Solute.calc_loggamma(
          T,
          p,
          X);
Real[nLfun - 1] molality=calc_mfromX(X);

algorithm

for i in 1:nLfun-1 loop
  ai[i] := 10^loggamma[i] * molality[i];
end for;

end calc_activity_i;
