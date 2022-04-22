within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Bromley.Reduced.Additional.Solvent.tautau;
function calc_loggammatautau "Helper function to calculate Gibbs derivative"

  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nLfun] X;

  output Real loggammatautau;

protected
  Real[nLfun - 1,nLfun - 1] W_ij=Bromley.Solvent.calc_W_ij(X);
  Real[nLfun-1,nLfun-1] log_a_ijtautau = calc_log_a_ijtautau(T,p,X);

algorithm

  loggammatautau :=sum(W_ij .* log_a_ijtautau);

end calc_loggammatautau;
