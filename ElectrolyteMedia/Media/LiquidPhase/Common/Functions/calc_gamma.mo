within ElectrolyteMedia.Media.LiquidPhase.Common.Functions;
function calc_gamma
  "Calculates activity coefficient of solutes and solvent at T and p"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nLfun] X;
  output SI.SpecificEnergy[nLfun] gamma;
algorithm
  if LiquidModelfun == Media.Common.Types.LiquidModel.Ideal then
    gamma :=calc_gamma_ideal();
  elseif LiquidModelfun == Media.Common.Types.LiquidModel.Debye then
    gamma :=GibbsExcess.Debye.calc_gamma(T,p,X);
  elseif LiquidModelfun == Media.Common.Types.LiquidModel.DebyeHueckel then
    gamma :=GibbsExcess.DebyeHueckel.calc_gamma(T,p,X);
  elseif LiquidModelfun == Media.Common.Types.LiquidModel.Bromley then
    gamma :=GibbsExcess.Bromley.calc_gamma(T,p,X);
  elseif LiquidModelfun == Media.Common.Types.LiquidModel.Pitzer then
    gamma :=GibbsExcess.Pitzer.calc_gamma(T,p,X);
  end if;

  annotation(smoothOrder=5);
end calc_gamma;
