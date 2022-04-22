within ElectrolyteMedia.Media.GasLiquidPhase.Common.Functions;
function calc_MMX
  output SI.MolarMass[nGfunfun+nLfunfun] MMX;

protected
  SI.MolarMass[nGfunfun] MMXg = GasFunctions.calc_MMX();
  SI.MolarMass[nLfunfun] MMXl = LiquidFunctions.calc_MMX();
algorithm
  MMX[1:nGfunfun] :=MMXg;// :=cat(1,MMXg,MMXl);
  MMX[1+nGfunfun:nGfunfun+nLfunfun] :=MMXl;
annotation(smoothOrder=5);
end calc_MMX;
