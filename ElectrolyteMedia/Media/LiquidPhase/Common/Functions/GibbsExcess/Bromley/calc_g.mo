within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Bromley;
function calc_g "Calculates specific Gibbs free energy with Bromley model"

  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nLfun] X;

  output SI.SpecificGibbsFreeEnergy[nLfun] g;

protected
  Real[nLfun] ai=calc_ai(T,p,X);
  SI.SpecificEntropy[nLfun] R = calc_RX();
algorithm

  for i in 1:nLfun-1 loop
    if ai[i] > 0 then
      g[i] := T *R[i]* log(ai[i]);
    end if;
  end for;

  if ai[nLfun] > 0 then
    g[nLfun] := T * R[nLfun]*log(ai[nLfun]);
  end if;

end calc_g;
