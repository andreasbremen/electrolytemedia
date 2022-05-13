within ElectrolyteMedia.Media.SolidLiquidPhase.Common.Functions;
function calc_MMX
  output SI.MolarMass[nSfunfun+nLfunfun] MMX;

protected
  SI.MolarMass[nSfunfun] MMXs = SolidFunctions.calc_MMX();
  SI.MolarMass[nLfunfun] MMXl = LiquidFunctions.calc_MMX();
algorithm
  MMX[1:nSfunfun] :=MMXs;// :=cat(1,MMXs,MMXl);
  MMX[1+nSfunfun:nSfunfun+nLfunfun] :=MMXl;

end calc_MMX;
