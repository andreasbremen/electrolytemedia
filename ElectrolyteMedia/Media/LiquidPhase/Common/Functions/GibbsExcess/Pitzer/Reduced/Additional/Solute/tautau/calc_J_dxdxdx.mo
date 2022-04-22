within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solute.tautau;
function calc_J_dxdxdx "Helper function to calculate Gibbs derivative"
  input Real x;
  output Real J_dxdxdx;
protected
  Real[23] d;
  Real z_dxdxdx;
algorithm
  z_dxdxdx := calc_z_dxdxdx(x);
  d := Pitzer.Solute.calc_d(x);

  J_dxdxdx :=0.5* z_dxdxdx *(d[1] - d[3]);

end calc_J_dxdxdx;
