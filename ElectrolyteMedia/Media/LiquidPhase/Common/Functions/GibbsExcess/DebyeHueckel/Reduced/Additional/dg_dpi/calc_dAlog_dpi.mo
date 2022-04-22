within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel.Reduced.Additional.dg_dpi;
function calc_dAlog_dpi "Helper function to calculate Gibbs derivative"
  input SI.Temperature T;
  input SI.Pressure p;
  output Real dAlog_dpi;

algorithm

    dAlog_dpi :=3/log(10)*calc_dAln_dpi(T, p);

end calc_dAlog_dpi;
