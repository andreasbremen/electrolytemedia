within ElectrolyteMedia.Media.LiquidPhase.Common.Functions;
function calc_v
  "Calculates specific volume of solute and solvent mixture at T and p"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[:] X;
  output SI.SpecificVolume v;

protected
  SI.SpecificVolume[nLfun] vi;
  SI.SpecificVolume[nLfun] vex;
algorithm
  vi[1:nLfun-1] :=Solutes.calc_v_i(T, p);
  vi[nLfun] :=IF97_R1_Tp.calc_v(T, p);

  if LiquidModelfun == Media.Common.Types.LiquidModel.DebyeHueckel then
    vex :=GibbsExcess.DebyeHueckel.calc_v(T,p,X);
  elseif LiquidModelfun == Media.Common.Types.LiquidModel.Bromley then
    vex :=GibbsExcess.Bromley.calc_v(T,p,X);
  else
    vex :=zeros(nLfun);
  end if;

  v:=X*(vi+vex);

  annotation(smoothOrder=5);
end calc_v;
