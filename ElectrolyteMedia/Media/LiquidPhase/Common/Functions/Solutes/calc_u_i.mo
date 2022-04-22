within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.Solutes;
function calc_u_i
  "Calculates specific internal energy of solutes at inifinite dilution based on HKF model"
  input SI.Temperature T;
  input SI.Pressure p;
  output SI.SpecificInternalEnergy[nLifun] u;
protected
  GibbsDerivs[nLifun] g = Reduced.GibbsDerivs(T,p);
  ThermoProperties_pT[nLifun] pro;
algorithm
  for i in 1:nLifun loop
    if datafun[i].MM > 0.001008 then
      pro[i] :=Modelica.Media.Common.ThermoFluidSpecial.gibbsToProps_pT(g[i]);
      u[i] :=pro[i].u;
    else
      u[i] :=0;
    end if;
  end for;
//   u :=pro.u;

  annotation(smoothOrder=5);
end calc_u_i;
