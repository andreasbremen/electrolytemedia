within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.Solutes;
function calc_cv_i
  "Calculates specific heat capacity cv of solutes at inifinite dilution based on HKF model"
  input SI.Temperature T;
  input SI.Pressure p;
  output SI.SpecificHeatCapacityAtConstantVolume [nLifun] cv;
protected
  GibbsDerivs[nLifun] g = Reduced.GibbsDerivs(T,p);
  ThermoProperties_pT[nLifun] pro;
algorithm
  for i in 1:nLifun loop
    if datafun[i].MM > 0.001008 then
      pro[i] :=Modelica.Media.Common.ThermoFluidSpecial.gibbsToProps_pT(g[i]);
      cv[i] :=pro[i].cv;
    else
      cv[i] :=0;
    end if;
  end for;

end calc_cv_i;
