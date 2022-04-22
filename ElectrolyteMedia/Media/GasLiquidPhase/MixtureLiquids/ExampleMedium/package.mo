within ElectrolyteMedia.Media.GasLiquidPhase.MixtureLiquids;
package ExampleMedium
extends Common.MixtureLiquid(userInterface(
    GasModel=ElectrolyteMedia.Media.Common.Types.GasModel.Ideal,
    datag=GasPhase.Common.SingleGasesData.H2O_CO2_O2,
    interactionl=LiquidPhase.Common.MixtureSolutesData.CarbonDioxide_NaCl,
    nR_GL=2,
    nu_GL=[-1,0,0,0,0,0,0,0,0,0,1; 0,1,0,-1,0,0,0,0,0,0,0],
    nR_L=3,
    nu_L=[-1,1,0,0,0,0,1,-1; 0,-1,1,0,0,0,1,0; 0,0,0,0,0,1,1,-1],
    mediumName="TestMedium",
      Tstart=298.15,
      pstart=100000,
    nG=3,
    gasSubstanceNames={"H2O,g","CO2,g","O2,g"},
    interactiong=GasPhase.Common.MixtureGasesData.H2O_CO2_O2,
    nL=8,
    LiquidModel=Media.Common.Types.LiquidModel.Ideal,
    liquidSubstanceNames={"CO2,l","HCO3-","CO3-2","Na+","Cl-","OH-","H+","H2O"},
    datal=LiquidPhase.Common.SolutesData.CarbonDioxide_NaCl,
    useLiquidMassFraction=false,
    useLiquidMolality=true,
    refml={0,0,0,1,1,0,0},
    useGasMassFraction=true,
      refXg={0,0.8,0.2},
    usePhaseMassFraction=true,
    usePhaseMoleFraction=false,
      refLiquidMassFraction=0.8));

end ExampleMedium;
