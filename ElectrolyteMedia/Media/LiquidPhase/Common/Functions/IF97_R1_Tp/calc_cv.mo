within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.IF97_R1_Tp;
function calc_cv
  "Calculates heat capacity cv of water with IF97 model and standard state reference"
  input SI.Temperature T;
  input SI.Pressure p;
  output SI.SpecificHeatCapacityAtConstantVolume cv;
protected
  GibbsDerivs g = Reduced.GibbsDerivs(T,p);
  ThermoProperties_pT pro;
algorithm
  pro :=Modelica.Media.Common.ThermoFluidSpecial.gibbsToProps_pT(g);
  cv :=pro.cv;

end calc_cv;
