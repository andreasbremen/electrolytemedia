within ElectrolyteMedia.Media.GasLiquidPhase.SaturatedLiquids;
package ExampleMedium
extends Common.SaturatedLiquid(userInterface(
    GasModel=ElectrolyteMedia.Media.Common.Types.GasModel.Ideal,
    datag=GasPhase.Common.SingleGasesData.H2O_array,
    interactionG=GasPhase.Common.MixtureGasesData.SingleGas,
    interactionL=LiquidPhase.Common.MixtureSolutesData.NaOH_dis_,
    nR_GL=1,
    nu_GL=[-1,0,0,0,0,1],
    nR_L=2,
    nu_L=[0,-1,0,-1,1; -1,-1,1,0,0],
    mediumName="TestMedium",
    Tstart=298.15,
    pstart=1000000,
    nG=1,
    gasSubstanceNames={"H2O,g"},
    nL=5,
    LiquidModel=ElectrolyteMedia.Media.Common.Types.LiquidModel.Pitzer,
    liquidSubstanceNames={"Na+","OH-","NaOH,aq","H+","H2O"},
    datal=LiquidPhase.Common.SolutesData.NaOH_dis_,
    useLiquidMassFraction=false,
    useLiquidMolality=true,
    refml={1,1,0,0},
    useGasMassFraction=true,
    refXg={1},
    usePhaseMassFraction=true,
    usePhaseMoleFraction=false));
end ExampleMedium;
