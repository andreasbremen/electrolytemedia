within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solute.tau;
function calc_J_dxdx "Helper function to calculate Gibbs derivative"
  input Real x;
  output Real J_dxdx;
protected
  Real[23] d;
  Real z_dxdx;
algorithm
  z_dxdx := calc_z_dxdx(x);
  d := Pitzer.Solute.calc_d(x);

  J_dxdx :=0.5* z_dxdx *(d[1] - d[3]);

end calc_J_dxdx;
