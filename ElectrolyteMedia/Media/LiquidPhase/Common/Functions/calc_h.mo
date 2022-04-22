within ElectrolyteMedia.Media.LiquidPhase.Common.Functions;
function calc_h
  "Calculates specific enthalpy of solute and solvent mixture at T and p"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nLfun] X;
  output SI.SpecificEnthalpy h;

protected
  SI.SpecificEnthalpy[nLfun] hi;
  SI.SpecificEnthalpy[nLfun] hex;
algorithm
  hi[1:nLfun-1] :=Solutes.calc_h_i(T, p);
  hi[nLfun] :=IF97_R1_Tp.calc_h(T, p);

  if LiquidModelfun == Media.Common.Types.LiquidModel.DebyeHueckel then
    hex :=GibbsExcess.DebyeHueckel.calc_h(T,p,X);
  elseif LiquidModelfun == Media.Common.Types.LiquidModel.Bromley then
    hex :=GibbsExcess.Bromley.calc_h(T,p,X);
  else
    hex :=zeros(nLfun);
  end if;

  h:=X*(hi+hex);

  annotation(smoothOrder=5);
end calc_h;
