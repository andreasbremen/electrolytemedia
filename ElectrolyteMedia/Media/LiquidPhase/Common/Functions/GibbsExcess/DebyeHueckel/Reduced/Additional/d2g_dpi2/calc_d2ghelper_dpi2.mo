within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel.Reduced.Additional.d2g_dpi2;
function calc_d2ghelper_dpi2 "Helper function to calculate Gibbs derivative"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction X[nLfun];
  output Real d2g_dpi2[nLfun-1];

protected
  Real I = calc_I(X);
algorithm

  for i in 1:nLfun-1 loop
    d2g_dpi2[i] :=1e-10*3.72*sqrt(I)*calc_d2Blog_dpi2(T, p);//datafun[i].a_0
  end for;

end calc_d2ghelper_dpi2;
