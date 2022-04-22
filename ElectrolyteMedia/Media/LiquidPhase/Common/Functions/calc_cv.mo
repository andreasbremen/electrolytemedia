within ElectrolyteMedia.Media.LiquidPhase.Common.Functions;
function calc_cv
  "Calculates specific heat capacity cv of solute and solvent mixture at T and p"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[:] X;
  output SI.SpecificHeatCapacityAtConstantVolume cv;

protected
  SI.SpecificHeatCapacityAtConstantVolume[nLfun] cvi;
  SI.SpecificHeatCapacityAtConstantVolume[nLfun] cvex;
algorithm
  cvi[1:nLfun-1] :=Solutes.calc_cv_i(T, p);
  cvi[nLfun] :=IF97_R1_Tp.calc_cv(T, p);

  if LiquidModelfun == Media.Common.Types.LiquidModel.DebyeHueckel then
    cvex :=GibbsExcess.DebyeHueckel.calc_cv(T,p,X);
  elseif LiquidModelfun == Media.Common.Types.LiquidModel.Bromley then
    cvex :=GibbsExcess.Bromley.calc_cv(T,p,X);
  else
    cvex :=zeros(nLfun);
  end if;

  cv :=X*(cvi+cvex);

end calc_cv;
