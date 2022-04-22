within ElectrolyteMedia.Media.LiquidPhase.Common.Functions;
function calc_g
  "Calculates specific Gibbs free energy of solute and solvent mixture at T and p"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[:] X;
  output SI.SpecificGibbsFreeEnergy g;

protected
  SI.SpecificGibbsFreeEnergy[nLfun] gi;
  SI.SpecificGibbsFreeEnergy[nLfun] gex;
  SI.SpecificGibbsFreeEnergy g_sym;
  SI.SpecificGibbsFreeEnergy gim;
  Real [nLfun] mi;
  Real [nLfun] Ri;
  Real [nLfun] Y;
  SI.SpecificEntropy R;
algorithm
  gi[1:nLfun-1] :=Solutes.calc_g_i(T, p);
  gi[nLfun] :=IF97_R1_Tp.calc_g(T, p);

  Y := calc_Y(X);
  mi[nLfun] :=Y[nLfun];
  mi[1:nLfun-1]:=calc_mfromX(X);
  Ri := calc_RX();

  if LiquidModelfun == Media.Common.Types.LiquidModel.Debye then
    gex :=GibbsExcess.Debye.calc_g(T,p,X);
  elseif LiquidModelfun == Media.Common.Types.LiquidModel.DebyeHueckel then
    gex :=GibbsExcess.DebyeHueckel.calc_g(T,p,X);
  elseif LiquidModelfun == Media.Common.Types.LiquidModel.Bromley then
    gex :=GibbsExcess.Bromley.calc_g(T,p,X);
  elseif LiquidModelfun == Media.Common.Types.LiquidModel.Pitzer then
    gex :=GibbsExcess.Pitzer.calc_g(T,p,X);
  else
    for i in 1:nLfun loop
      if mi[i] > 0 then
        gex[i] :=T*Ri[i]*log(mi[i]);//zeros(nLfun);//
      else
        gex[i] :=0;
      end if;
    end for;
  end if;

  g :=X*gi+ X*gex;

end calc_g;
