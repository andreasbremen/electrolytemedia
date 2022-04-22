within ElectrolyteMedia.Media.LiquidPhase.Common;
package Functions "Functions to calculate thermodynamic properties for liquid phase"

   constant Integer nLfun=2 "Number of species";
   final constant Integer nLifun(min= 1) = nLfun-1 "Number of species -1";
   constant DataRecord[nLfun-1] datafun "DataRecord containing thermodynamic properties";
   constant LiquidInteractionDataRecord interactionfun "DataRecord containing thermodynamic properties for interaction between species";
   constant Media.Common.Types.LiquidModel LiquidModelfun "Type of liquid model used for calculations";















end Functions;
