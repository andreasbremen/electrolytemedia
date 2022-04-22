within ElectrolyteMedia.Media.GasLiquidPhase.Common;
package Functions

  constant Integer nGfunfun(min=1)=1;
  constant DataRecordG[nGfunfun] dataGfunfun;
constant GasInteractionDataRecord interactionGfunfun;
  constant Media.Common.Types.GasModel GasModelfunfun;

  constant Integer nLfunfun=2;
  final constant Integer nLifunfun = nLfunfun-1;
  constant Common.DataRecordL[nLfunfun - 1] dataLfunfun;
  constant Common.LiquidInteractionDataRecord interactionLfunfun;
  constant Media.Common.Types.LiquidModel LiquidModelfunfun;

end Functions;
