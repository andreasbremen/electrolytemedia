within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel.Reduced.Additional.d2g_dtau2;
function calc_d2ai_dtau2 "Helper function to calculate Gibbs derivative"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nLfun] X;
  output Real[nLfun-1] d2ai_dtau2;

protected
  Real mol[nLfun-1] = calc_mfromX(X);
  Real loggamma[nLfun - 1]=calc_loggamma_i(
                  T,
                  p,
                  X);
  Real dloggamma_dtau[nLfun-1] = dg_dtau.calc_dloggamma_dtau(T,p,X);
  Real d2loggamma_dtau2[nLfun-1] = d2g_dtau2.calc_d2loggamma_dtau2(T,p,X);

algorithm

  for i in 1:nLfun-1 loop
    d2ai_dtau2[i] := mol[i]*log(10)*(log(10)*10^loggamma[i]*(dloggamma_dtau[i])^2+10^loggamma[i]*d2loggamma_dtau2[i]);
  end for;

end calc_d2ai_dtau2;
