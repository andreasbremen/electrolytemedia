within ElectrolyteMedia.Media.GasPhase.MixtureGases;
package CO2_H2O_IG "Ideal Gas EOS: CO2, H2O"
  extends Common.MixtureGas(
    userInterface(
      nG=2,
      GasModel=ElectrolyteMedia.Media.Common.Types.GasModel.Ideal,
      gasSubstanceNames={"CO2","H2O"},
      datag=Common.SingleGasesData.CO2_H2O,
      fluidconstants=Common.FluidData.CO2_H2O,
      mediumName="Ideal gas with CO2 and H2O",
      Tstart=298.15,
      pstart=100000,
      useMassFraction=true,
      refXg={0.6,0.4}));
end CO2_H2O_IG;
