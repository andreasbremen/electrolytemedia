within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.Solutes;
function calc_s_i
  "Calculates specific entropy of solutes at inifinite dilution based on HKF model"
  input SI.Temperature T;
  input SI.Pressure p;
  output SI.SpecificEntropy[nLifun] s;
protected
  GibbsDerivs[nLifun] g = Reduced.GibbsDerivs_(T,p);
  ThermoProperties_pT[nLifun] pro;
algorithm
  for i in 1:nLifun loop
    if datafun[i].MM > 0.001008 then
      pro[i] :=Modelica.Media.Common.ThermoFluidSpecial.gibbsToProps_pT(g[i]);
      s[i] :=pro[i].s;
    else
      s[i] :=0;
    end if;
//     s[i] :=pro[i].s;
  end for;

end calc_s_i;
