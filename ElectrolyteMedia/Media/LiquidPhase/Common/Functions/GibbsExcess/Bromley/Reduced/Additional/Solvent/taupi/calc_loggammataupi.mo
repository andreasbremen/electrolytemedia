within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Bromley.Reduced.Additional.Solvent.taupi;
function calc_loggammataupi "Helper function to calculate Gibbs derivative"

  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nLfun] X;

  output Real loggammataupi;

protected
  Real[nLfun - 1,nLfun - 1] W_ij=Bromley.Solvent.calc_W_ij(X);
  Real[nLfun-1,nLfun-1] log_a_ijtaupi = calc_log_a_ijtaupi(T,p,X);

algorithm

  loggammataupi :=sum(W_ij .* log_a_ijtaupi);

end calc_loggammataupi;
