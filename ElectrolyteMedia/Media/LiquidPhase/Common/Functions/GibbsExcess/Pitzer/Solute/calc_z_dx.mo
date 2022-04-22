within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Solute;
function calc_z_dx "Unsymmetrical mixing, Pitzer p.124"
  input Real x;
  output Real der_z;
algorithm
  if x < 1 then
    der_z :=4/5*x^(-4/5);
  else
    der_z :=-40/90*x^(-11/10);
  end if;
end calc_z_dx;
