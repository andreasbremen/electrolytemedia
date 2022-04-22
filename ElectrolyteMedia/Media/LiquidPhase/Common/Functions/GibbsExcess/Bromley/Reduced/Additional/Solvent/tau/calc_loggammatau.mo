within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Bromley.Reduced.Additional.Solvent.tau;
function calc_loggammatau "Helper function to calculate Gibbs derivative"

  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nLfun] X;

  output Real loggammatau;

protected
  Real[nLfun - 1,nLfun - 1] W_ij=Bromley.Solvent.calc_W_ij(X);
  Real[nLfun-1,nLfun-1] log_a_ijtau = calc_log_a_ijtau(T,p,X);

algorithm

  loggammatau :=sum(W_ij .* log_a_ijtau);

end calc_loggammatau;
