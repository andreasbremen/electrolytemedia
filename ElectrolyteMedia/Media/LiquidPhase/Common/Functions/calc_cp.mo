within ElectrolyteMedia.Media.LiquidPhase.Common.Functions;
function calc_cp
  "Calculates specific heat capacity cp of solute and solvent mixture at T and p"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[:] X;
  output SI.SpecificHeatCapacityAtConstantPressure cp;

protected
  SI.SpecificHeatCapacityAtConstantPressure[nLfun] cpi;
  SI.SpecificHeatCapacityAtConstantPressure[nLfun] cpex;
algorithm
  cpi[1:nLfun-1] :=Solutes.calc_cp_i(T, p);
  cpi[nLfun] :=IF97_R1_Tp.calc_cp(T, p);

  if LiquidModelfun == Media.Common.Types.LiquidModel.DebyeHueckel then
    cpex :=GibbsExcess.DebyeHueckel.calc_cp(T,p,X);
  elseif LiquidModelfun == Media.Common.Types.LiquidModel.Bromley then
    cpex :=GibbsExcess.Bromley.calc_cp(T,p,X);
  else
    cpex :=zeros(nLfun);
  end if;

  cp :=X*(cpi+cpex);

end calc_cp;
