within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel.Reduced.Additional.dg_dtau;
function calc_dAlog_dtau "Helper function to calculate Gibbs derivative"
  input SI.Temperature T;
  input SI.Pressure p;
  output Real dAlog_dT;

algorithm

    dAlog_dT :=3/log(10)*calc_dAln_dtau(T, p);

end calc_dAlog_dtau;
