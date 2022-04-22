within ElectrolyteMedia.Media.GasLiquidPhase.SaturatedLiquids;
package NaCl
extends Common.SaturatedLiquid(
     userInterface(
    GasModel=ElectrolyteMedia.Media.Common.Types.GasModel.Ideal,
    datag=GasPhase.Common.SingleGasesData.H2O_array,
      nR_GL=1,
      nu_GL=[1,0,0,-1],
      nR_L=0,
    interactionG=GasPhase.Common.MixtureGasesData.SingleGas,
      interactionL=LiquidPhase.Common.MixtureSolutesData.NaCl,
       mediumName="TestMedium",
      Tstart=298.15,
      pstart=1000000,
      nG=1,
      gasSubstanceNames={"H2O,g"},
      nL=3,
    LiquidModel=ElectrolyteMedia.Media.Common.Types.LiquidModel.Pitzer,
      liquidSubstanceNames={"Na+","Cl-","H2O"},
      datal=LiquidPhase.Common.SolutesData.NaCl,
    useLiquidMassFraction=false,
      useLiquidMolality=true,
      refml={1,1},
      useGasMassFraction=true,
      refXg={1},
      usePhaseMassFraction=true,
    usePhaseMoleFraction=false));

end NaCl;
