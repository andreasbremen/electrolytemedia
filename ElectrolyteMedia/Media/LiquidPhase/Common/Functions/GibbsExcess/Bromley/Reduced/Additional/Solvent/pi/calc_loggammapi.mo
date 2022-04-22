within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Bromley.Reduced.Additional.Solvent.pi;
function calc_loggammapi "Helper function to calculate Gibbs derivative"

  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nLfun] X;

  output Real loggammapi;

protected
  Real[nLfun - 1,nLfun - 1] W_ij=Bromley.Solvent.calc_W_ij(X);
  Real[nLfun-1,nLfun-1] log_a_ijpi = calc_log_a_ijpi(T,p,X);

algorithm

  loggammapi :=sum(W_ij .* log_a_ijpi);

end calc_loggammapi;
