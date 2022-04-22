within ElectrolyteMedia.Media.GasPhase.Common.MixtureGas;
package Functions "Extends Common.Functions and passes data"
  extends Common.Functions(
    nGfun = nX,
    interactionfun = interaction,
    GasModelfun=GasModel,
    datafun = data);

end Functions;
