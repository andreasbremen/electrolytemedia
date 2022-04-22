within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel.Reduced.Additional.d2g_dpi2;
function calc_dhhelper_dpi "Helper function to calculate Gibbs derivative"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction X[nLfun];
  output Real dh_dpi[nLfun - 1];

protected
  Real[nLfun-1] d2fhelper_dpi2 = calc_d2fhelper_dpi2(T,p,X);
  Real[nLfun-1] ghelper = dg_dpi.calc_ghelper(T,p,X);
  Real[nLfun-1] d2ghelper_dpi2= calc_d2ghelper_dpi2(T,p,X);
  Real[nLfun-1] fhelper= dg_dpi.calc_fhelper(T,p,X);

algorithm

  for i in 1:nLfun-1 loop
    dh_dpi[i] := d2fhelper_dpi2[i]*ghelper[i]-(d2ghelper_dpi2[i]*fhelper[i]);
  end for;

end calc_dhhelper_dpi;
