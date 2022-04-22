within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.IF97_R1_Tp;
function calc_f
  "Calculates specific Helmholtz energy of water with IF97 model and standard state reference"
  input SI.Temperature T;
  input SI.Pressure p;
  output SI.SpecificGibbsFreeEnergy f;
protected
  GibbsDerivs g = Reduced.GibbsDerivs(T,p);
  ThermoProperties_pT pro;
algorithm
  pro :=Modelica.Media.Common.ThermoFluidSpecial.gibbsToProps_pT(g);
  f :=StdRefH2O.a_tr + pro.u - (T * (pro.s + StdRefH2O.s_tr) - StdRefH2O.T_tr * StdRefH2O.s_tr);

end calc_f;
