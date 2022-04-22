within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel.Reduced;
function calc_g
  "Calculates excess reduced Gibbs free energy at infinite dilution of aqueous species"

  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nLfun] X;
  output Real[nLfun] g;

protected
  Real[nLfun] a=calc_activity(T,p,X);
algorithm
  for i in 1:nLfun loop
    if a[i] > 0 then
      g[i] := log(a[i]);
    end if;
  end for;

end calc_g;
