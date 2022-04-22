within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel.Reduced.Additional.d2g_dtaudpi;
function calc_dkhelper_dpi "Helper function to calculate Gibbs derivative"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction X[nLfun];
  output Real dk_dpi[nLfun - 1];

protected
  Real[nLfun-1] ghelper= dg_dpi.calc_ghelper(T,p,X);
  Real[nLfun-1] dghelper_dpi=dg_dpi.calc_dghelper_dpi(T,p,X);

algorithm

  for i in 1:nLfun-1 loop
    dk_dpi[i] := 2*ghelper[i]*dghelper_dpi[i];
  end for;

end calc_dkhelper_dpi;
