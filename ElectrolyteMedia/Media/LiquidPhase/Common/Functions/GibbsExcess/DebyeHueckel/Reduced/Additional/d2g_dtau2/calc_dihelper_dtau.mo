within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel.Reduced.Additional.d2g_dtau2;
function calc_dihelper_dtau "Helper function to calculate Gibbs derivative"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction X[nLfun];
  output Real di_dT[nLfun - 1];

protected
  Real[nLfun-1] ghelper= dg_dtau.calc_ghelper(T,p,X);
  Real[nLfun-1] dghelper_dtau=dg_dtau.calc_dghelper_dtau(T,p,X);

algorithm

  for i in 1:nLfun-1 loop
    di_dT[i] := 2*ghelper[i]*dghelper_dtau[i];
  end for;

end calc_dihelper_dtau;
