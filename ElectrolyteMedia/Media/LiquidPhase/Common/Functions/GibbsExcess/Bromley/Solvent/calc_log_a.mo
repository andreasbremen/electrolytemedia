within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Bromley.Solvent;
function calc_log_a
  "calculates decadic log of activity of water based on Bromleys method with extension to multielectrolyte systems by Meissner and Kusik"

  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nLfun] X;

  output Real log_a;
protected
  Real W_ij[nLfun - 1,nLfun - 1]=calc_W_ij(X);
  Real r=Solvent.calc_r(X);
  Real[nLfun - 1,nLfun - 1] log_a_ij=calc_log_a_ij(T,p,X);
algorithm

  log_a :=sum(W_ij .* log_a_ij) + r;
annotation(smoothOrder=20);
end calc_log_a;
