within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel.Reduced.Additional.dg_dpi;
function calc_dai_dpi "Helper function to calculate Gibbs derivative"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction X[nLfun];
  output Real dai_dpi[nLfun-1];

protected
  Real[nLfun-1] m_i = calc_mfromX(X);
  Real[nLfun - 1] loggamma=calc_loggamma_i(
                  T,
                  p,
                  X);
  Real[nLfun - 1] dloggamma_dpi=calc_dloggamma_dpi(T,p,X);

algorithm
  for i in 1:nLfun-1 loop
  dai_dpi[i] := m_i[i] * log(10) * 10^(loggamma[i]) * dloggamma_dpi[i];
  end for;

end calc_dai_dpi;
