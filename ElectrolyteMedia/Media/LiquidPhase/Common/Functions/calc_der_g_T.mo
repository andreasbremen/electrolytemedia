within ElectrolyteMedia.Media.LiquidPhase.Common.Functions;
function calc_der_g_T
  "Calculates temperature derivative of specific Gibbs free energy of solute and solvent mixture at T and p"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[:] X;
  output Real gdT(unit="J/(kg.K)");

protected
  Real[nLfun] gidT(unit="J/(kg.K)");
  Real[nLfun] gexdT(unit="J/(kg.K)");
algorithm
  gidT[1:nLfun-1] :=Solutes.calc_der_g_T(T, p);
  gidT[nLfun] :=IF97_R1_Tp.calc_der_g_T(T, p);

  if LiquidModelfun == Media.Common.Types.LiquidModel.DebyeHueckel then
    gexdT :=GibbsExcess.DebyeHueckel.calc_der_g_T(
      T,
      p,
      X);
  elseif LiquidModelfun == Media.Common.Types.LiquidModel.Bromley then
    gexdT :=GibbsExcess.Bromley.calc_der_g_T(
      T,
      p,
      X);
  else
    gexdT :=zeros(nLfun);
  end if;

  gdT :=X*(gidT+gexdT);

end calc_der_g_T;
