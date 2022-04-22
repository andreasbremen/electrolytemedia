within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel.Solvent;
function calc_sigma
  input Real x;
  output Real sigma;
algorithm
  if x > 0 then
    sigma :=3/(x)^3*((1 + x) - 2*log(1 + x) - 1/(1 + x));
  else
    sigma :=1;
  end if;
end calc_sigma;
