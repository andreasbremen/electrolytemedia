within ElectrolyteMedia.Media.LiquidPhase.Common.Functions;
function calc_u
  "Calculates specific internal energy of solute and solvent mixture at T and p"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nLfun] X;
  output SI.SpecificInternalEnergy u;

protected
  SI.SpecificInternalEnergy[nLfun] ui;
  SI.SpecificInternalEnergy[nLfun] uex;
algorithm
  ui[1:nLfun-1] :=Solutes.calc_u_i(T, p);
  ui[nLfun] :=IF97_R1_Tp.calc_u(T, p);

  if LiquidModelfun == Media.Common.Types.LiquidModel.DebyeHueckel then
    uex :=GibbsExcess.DebyeHueckel.calc_u(T,p,X);
  elseif LiquidModelfun == Media.Common.Types.LiquidModel.Bromley then
    uex :=GibbsExcess.Bromley.calc_u(T,p,X);
  else
    uex :=zeros(nLfun);
  end if;

  u :=X*(ui+uex);

  annotation(smoothOrder=5);
end calc_u;
