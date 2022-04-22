within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel.Reduced.Additional.d2g_dtaudpi;
function calc_d2Alog_dtaudpi "Helper function to calculate Gibbs derivative"
  input SI.Temperature T;
  input SI.Pressure p;
  output Real d2Alog_dtaudpi;

algorithm

  d2Alog_dtaudpi :=3/log(10)*calc_d2Aln_dtaudpi(T, p);

end calc_d2Alog_dtaudpi;
