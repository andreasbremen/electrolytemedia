within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.Solutes;
function calc_v_i
  "Calculates specific volume of solutes at inifinite dilution based on HKF model"
  input SI.Temperature T;
  input SI.Pressure p;
  output SI.SpecificVolume[nLifun] v;
protected
  GibbsDerivs[nLifun] g = Reduced.GibbsDerivs(T,p);
  ThermoProperties_pT[nLifun] pro;
algorithm
  for i in 1:nLifun loop
    if datafun[i].MM > 0.001008 then
      pro[i] :=Modelica.Media.Common.ThermoFluidSpecial.gibbsToProps_pT(g[i]);
      v[i] :=1/pro[i].d;
    else
      v[i] :=0;
    end if;
  end for;
annotation(smoothOrder=2);
end calc_v_i;
