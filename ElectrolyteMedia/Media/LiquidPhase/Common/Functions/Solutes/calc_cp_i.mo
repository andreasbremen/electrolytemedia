within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.Solutes;
function calc_cp_i
  "Calculates specific heat capacity cp of solutes at inifinite dilution based on HKF model"
  input SI.Temperature T;
  input SI.Pressure p;
  output SI.SpecificHeatCapacityAtConstantPressure [nLifun] cp;
protected
  GibbsDerivs[nLifun] g = Reduced.GibbsDerivs(T,p);
  ThermoProperties_pT[nLifun] pro;
algorithm
  for i in 1:nLifun loop
    if datafun[i].MM > 0.001008 then
      pro[i] :=Modelica.Media.Common.ThermoFluidSpecial.gibbsToProps_pT(g[i]);
      cp[i] :=pro[i].cp;
    else
      cp[i] :=0;
    end if;
  end for;

end calc_cp_i;
