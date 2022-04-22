within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.IF97_R1_Tp;
function calc_h
  "Calculates specific enthalpy of water with IF97 model and standard state reference"
  input SI.Temperature T;
  input SI.Pressure p;
  output SI.SpecificEnthalpy h;
protected
  GibbsDerivs g = Reduced.GibbsDerivs(T,p);
  ThermoProperties_pT pro;
algorithm
  pro :=Modelica.Media.Common.ThermoFluidSpecial.gibbsToProps_pT(g);
  h :=pro.h + StdRefH2O.h_tr;

  annotation(smoothOrder=5);
end calc_h;
