within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel.Reduced.Additional.dg_dpi;
function calc_dghelper_dpi "Helper function to calculate Gibbs derivative"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction X[nLfun];
  output Real dg_dpi[nLfun-1];

protected
  Real I = calc_I(X);
algorithm

  for i in 1:nLfun-1 loop
    dg_dpi[i] :=1e-10*3.72*sqrt(I)*calc_dBlog_dpi(T, p);//datafun[i].a_0
  end for;

end calc_dghelper_dpi;
