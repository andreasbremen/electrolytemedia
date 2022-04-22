within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.IF97_R1_Tp;
function calc_v "Calculates specific volume from temperature and pressure"
  input SI.Temperature T;
  input SI.Pressure p;
  output SI.SpecificVolume v;
protected
  GibbsDerivs g = Reduced.GibbsDerivs(T,p);
  ThermoProperties_pT pro;
algorithm
  pro :=Modelica.Media.Common.ThermoFluidSpecial.gibbsToProps_pT(g);
  v :=1/pro.d;
  annotation(smoothOrder = 5);
end calc_v;
