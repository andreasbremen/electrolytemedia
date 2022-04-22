within ElectrolyteMedia.Media.GasPhase.MixtureGases;
package CO2_H2O_PR "Peng-Robinson EOS: CO2, H2O"
  extends Common.MixtureGas(
    userInterface(
      nG=2,
      GasModel=ElectrolyteMedia.Media.Common.Types.GasModel.PengRobinson,
      gasSubstanceNames={"CO2","H2O"},
      datag=Common.SingleGasesData.CO2_H2O,
    interactiong=Common.MixtureGasesData.CarbonationGas,
      fluidconstants=Common.FluidData.CO2_H2O,
      mediumName="Peng Robinson gas with CO2 and H2O",
    Tstart=298.15,
    pstart=100000,
      useMassFraction=true,
      refXg={0.6,0.4}));
end CO2_H2O_PR;
