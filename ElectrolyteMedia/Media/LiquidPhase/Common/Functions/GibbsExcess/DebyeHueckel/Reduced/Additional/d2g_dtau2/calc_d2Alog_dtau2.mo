within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel.Reduced.Additional.d2g_dtau2;
function calc_d2Alog_dtau2 "Helper function to calculate Gibbs derivative"
  input SI.Temperature T;
  input SI.Pressure p;
  output Real d2Alog_dT2;

algorithm

  d2Alog_dT2 :=3/log(10)*calc_d2Aln_dtau2(T, p);

end calc_d2Alog_dtau2;
