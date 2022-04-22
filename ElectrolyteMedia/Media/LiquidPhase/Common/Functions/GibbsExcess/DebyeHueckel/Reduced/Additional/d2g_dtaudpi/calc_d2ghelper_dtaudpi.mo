within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel.Reduced.Additional.d2g_dtaudpi;
function calc_d2ghelper_dtaudpi "Helper function to calculate Gibbs derivative"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction X[nLfun];
  output Real d2g_dtaudpi[nLfun-1];

protected
  Real I = calc_I(X);
algorithm

  for i in 1:nLfun-1 loop
    d2g_dtaudpi[i] :=1e-10*3.72*sqrt(I)*calc_d2Blog_dtaudpi(T, p);//datafun[i].a_0
  end for;

end calc_d2ghelper_dtaudpi;
