within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.IF97_R1_Tp;
function calc_der_h_T
  "Calculates temperature derivative of specific enthalpy of water with IF97 model and standard state reference"
  input SI.Temperature T;
  input SI.Pressure p;
  output Real h_dT(unit="J/(kg.K)");
protected
  GibbsDerivs g = Reduced.GibbsDerivs(T,p);
  ThermoProperties_pT pro;
algorithm
  pro :=Modelica.Media.Common.ThermoFluidSpecial.gibbsToProps_pT(g);

  h_dT :=pro.cp;
end calc_der_h_T;
