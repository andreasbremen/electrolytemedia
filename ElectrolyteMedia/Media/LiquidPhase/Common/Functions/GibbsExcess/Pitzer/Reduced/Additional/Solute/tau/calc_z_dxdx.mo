within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solute.tau;
function calc_z_dxdx "Helper function to calculate Gibbs derivative"
  input Real x;
  output Real der_z_dxdx;
algorithm
  if x < 1 then
    der_z_dxdx :=-16/25*x^(-9/5);
  else
    der_z_dxdx :=22/45*x^(-21/10);
  end if;
end calc_z_dxdx;
