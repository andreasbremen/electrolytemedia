within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel.Reduced.Additional.d2g_dtau2;
function calc_dhhelper_dtau "Helper function to calculate Gibbs derivative"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction X[nLfun];
  output Real dh_dtau[nLfun - 1];

protected
  Real[nLfun-1] d2fhelper_dtau2 = calc_d2fhelper_dtau2(T,p,X);
  Real[nLfun-1] ghelper = dg_dtau.calc_ghelper(T,p,X);
  Real[nLfun-1] d2ghelper_dtau2= calc_d2ghelper_dtau2(T,p,X);
  Real[nLfun-1] fhelper= dg_dtau.calc_fhelper(T,p,X);

algorithm

  for i in 1:nLfun-1 loop
    dh_dtau[i] := d2fhelper_dtau2[i]*ghelper[i]-(d2ghelper_dtau2[i]*fhelper[i]);
  end for;

end calc_dhhelper_dtau;
