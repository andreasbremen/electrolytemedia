within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.IF97_R1_Tp;
function calc_der_s_T
  "Calculates temperature derivative of specific entropy of water with IF97 model and standard state reference"
  input SI.Temperature T;
  input SI.Pressure p;
  output Real s_dT(unit="J/(kg.K2)");
protected
  GibbsDerivs g = Reduced.GibbsDerivs(T,p);
  ThermoProperties_pT pro;
algorithm
  pro :=Modelica.Media.Common.ThermoFluidSpecial.gibbsToProps_pT(g);

  s_dT :=pro.cp/T;
end calc_der_s_T;
