within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Bromley.Reduced.Additional.Solvent.pipi;
function calc_loggammapipi "Helper function to calculate Gibbs derivative"

  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nLfun] X;

  output Real loggammapipi;

protected
  Real[nLfun - 1,nLfun - 1] W_ij=Bromley.Solvent.calc_W_ij(X);
  Real[nLfun-1,nLfun-1] log_a_ijpipi = calc_log_a_ijpipi(T,p,X);

algorithm

  loggammapipi :=sum(W_ij .* log_a_ijpipi);

end calc_loggammapipi;
