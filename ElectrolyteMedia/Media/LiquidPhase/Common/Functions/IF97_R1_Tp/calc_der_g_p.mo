within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.IF97_R1_Tp;
function calc_der_g_p
  "Calculates pressure derivative of specific Gibbs free energy of water with IF97 model and standard state reference"
  input SI.Temperature T;
  input SI.Pressure p;
  output Real gibbs_dp(unit="J/(kg.Pa)");
protected
  GibbsDerivs g = Reduced.GibbsDerivs(T,p);
  ThermoProperties_pT pro;
algorithm
  pro :=Modelica.Media.Common.ThermoFluidSpecial.gibbsToProps_pT(g);
// pro.hdp - T * sdp
  gibbs_dp := 1/pro.d;
end calc_der_g_p;
