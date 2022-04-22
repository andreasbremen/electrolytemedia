within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel;
function calc_cp
  "Calculates specific excess heat capacity of liquid species with Debye Hückel model"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nLfun] X;
  output SI.SpecificHeatCapacityAtConstantPressure[nLfun] cp;
protected
  GibbsDerivs[nLfun] g = Reduced.GibbsDerivs(T,p,X);
  ThermoProperties_pT[nLfun] pro;
algorithm
  for i in 1:nLfun-1 loop
    if datafun[i].z <> 0 and X[i] > 0 then
      pro[i] :=Modelica.Media.Common.ThermoFluidSpecial.gibbsToProps_pT(g[i]);
      cp[i] :=pro[i].cp;
    end if;
  end for;

end calc_cp;
