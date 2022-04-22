within ElectrolyteMedia.Media.LiquidPhase.Common.Functions;
function calc_s
  "Calculates specific entropy of solute and solvent mixture at T and p"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nLfun] X;
  output SI.SpecificEntropy s;

protected
  SI.SpecificEntropy[nLfun] si;
  SI.SpecificEntropy[nLfun] sex;
 // SI.SpecificEntropy ssym;
algorithm
  si[1:nLfun-1] :=Solutes.calc_s_i(T, p);
  si[nLfun] :=IF97_R1_Tp.calc_s(T, p);
  //ssym :=Solutes.Reduced.Sym.calc_s_sym(T,p,X);

  if LiquidModelfun == Media.Common.Types.LiquidModel.DebyeHueckel then
    sex :=GibbsExcess.DebyeHueckel.calc_s(T,p,X);
  elseif LiquidModelfun == Media.Common.Types.LiquidModel.Bromley then
    sex :=GibbsExcess.Bromley.calc_s(T,p,X);
  else
    sex :=zeros(nLfun);
  end if;

  s :=X*(si+sex);//+ssym;

  annotation(smoothOrder=5);
end calc_s;
