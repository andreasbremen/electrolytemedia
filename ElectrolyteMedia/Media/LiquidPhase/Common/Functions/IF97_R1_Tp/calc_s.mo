within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.IF97_R1_Tp;
function calc_s
  "Calculates specific entropy of water with IF97 model and standard state reference"
  input SI.Temperature T;
  input SI.Pressure p;
  output SI.SpecificEntropy s;
protected
  GibbsDerivs g = Reduced.GibbsDerivs(T,p);
  ThermoProperties_pT pro;
algorithm
  pro :=Modelica.Media.Common.ThermoFluidSpecial.gibbsToProps_pT(g);
  s :=pro.s + StdRefH2O.s_tr;

end calc_s;
