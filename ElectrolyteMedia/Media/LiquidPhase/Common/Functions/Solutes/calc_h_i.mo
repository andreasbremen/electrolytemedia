within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.Solutes;
function calc_h_i
  "Calculates specific enthalpy of solutes at inifinite dilution based on HKF model"
  input SI.Temperature T;
  input SI.Pressure p;
  output SI.SpecificEnthalpy[nLifun] h;
protected
  GibbsDerivs[nLifun] g = Reduced.GibbsDerivs(T,p);
  ThermoProperties_pT[nLifun] pro;
algorithm
  for i in 1:nLifun loop
    if datafun[i].MM > 0.001008 then
      pro[i] :=Modelica.Media.Common.ThermoFluidSpecial.gibbsToProps_pT(g[i]);
      h[i] :=pro[i].h;
    else
      h[i] :=0;
    end if;
  end for;

  annotation(smoothOrder=5);
end calc_h_i;
