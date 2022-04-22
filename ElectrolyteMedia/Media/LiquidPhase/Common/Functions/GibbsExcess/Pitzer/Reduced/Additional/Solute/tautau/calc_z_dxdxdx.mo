within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solute.tautau;
function calc_z_dxdxdx "Helper function to calculate Gibbs derivative"
  input Real x;
  output Real der_z_tautau;
algorithm
  if x < 1 then
    der_z_tautau :=144/125*x^(-14/5);
  else
    der_z_tautau :=-77/75*x^(-31/10);
  end if;
end calc_z_dxdxdx;
