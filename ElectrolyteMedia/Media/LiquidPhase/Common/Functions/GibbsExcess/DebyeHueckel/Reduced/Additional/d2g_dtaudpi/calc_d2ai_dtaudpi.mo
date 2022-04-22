within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel.Reduced.Additional.d2g_dtaudpi;
function calc_d2ai_dtaudpi "Helper function to calculate Gibbs derivative"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction X[nLfun];
  output Real d2ai_dtaudpi[nLfun - 1];

protected
  Real[nLfun-1] m_i = calc_mfromX(X);
  Real[nLfun - 1] loggamma=calc_loggamma_i(
                  T,
                  p,
                  X);
  Real[nLfun - 1] dloggamma_dtau=dg_dtau.calc_dloggamma_dtau(T,p,X);
  Real[nLfun - 1] dloggamma_dpi=dg_dpi.calc_dloggamma_dpi(T,p,X);
  Real[nLfun-1] d2loggamma_dtaudpi = calc_d2loggamma_dtaudpi(T,p,X);

algorithm

  for i in 1:nLfun-1 loop
    d2ai_dtaudpi[i] :=m_i[i]*log(10)*(log(10)*10^loggamma[i]*dloggamma_dpi[i]*dloggamma_dtau[i]+10^loggamma[i]*d2loggamma_dtaudpi[i]);
  end for;

end calc_d2ai_dtaudpi;
