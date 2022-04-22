within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel.Reduced.Additional.d2g_dpi2;
function calc_d2ai_dpi2 "Helper function to calculate Gibbs derivative"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nLfun] X;
  output Real[nLfun-1] d2ai_dpi2;

protected
  Real mol[nLfun-1] = calc_mfromX(X);
  Real loggamma[nLfun - 1]=calc_loggamma_i(
                  T,
                  p,
                  X);
  Real dloggamma_dpi[nLfun-1] = dg_dpi.calc_dloggamma_dpi(T,p,X);
  Real d2loggamma_dpi2[nLfun-1] = d2g_dpi2.calc_d2loggamma_dpi2(T,p,X);

algorithm

  for i in 1:nLfun-1 loop
    d2ai_dpi2[i] := mol[i]*log(10)*(log(10)*10^loggamma[i]*(dloggamma_dpi[i])^2+10^loggamma[i]*d2loggamma_dpi2[i]);
  end for;

end calc_d2ai_dpi2;
