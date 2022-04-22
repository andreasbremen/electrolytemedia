within ElectrolyteMedia.Media.GasPhase.Common;
package Functions "Functions to calculate thermodynamic properties"

  constant Integer nGfun(min=1)=1 "Number of gases";
  //constant Common.GasDataRecordIG[nGfun] dataIGfun "DataRecord containing thermodynamic properties for ideal gases";
  constant GasDataRecord[nGfun] datafun
    "DataRecord containing thermodynamic properties for Peng-Robinson gases";
  constant Common.InteractionDataRecord interactionfun "DataRecord containing thermodynamic properties for interaction between species";
  constant Media.Common.Types.GasModel GasModelfun "Type of gas model used for calculations";



















end Functions;
