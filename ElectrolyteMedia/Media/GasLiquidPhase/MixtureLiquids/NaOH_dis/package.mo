within ElectrolyteMedia.Media.GasLiquidPhase.MixtureLiquids;
package NaOH_dis
extends Common.MixtureLiquid(userInterface(
    GasModel=ElectrolyteMedia.Media.Common.Types.GasModel.PengRobinson,
    datag=GasPhase.Common.SingleGasesData.H2O_array,
    interactiong=GasPhase.Common.MixtureGasesData.SingleGas,
    interactionl=LiquidPhase.Common.MixtureSolutesData.NaOH_dis,
    nR_GL=1,
    nu_GL=[-1,0,0,0,1],
    nR_L=1,
    nu_L=[1,1,-1,0],
    mediumName="TestMedium",
    Tstart=298.15,
    pstart=1000000,
    nG=1,
    gasSubstanceNames={"H2O,g"},
    nL=4,
    LiquidModel=ElectrolyteMedia.Media.Common.Types.LiquidModel.Pitzer,
    liquidSubstanceNames={"Na+","OH-","NaOH,aq","H2O"},
    datal=LiquidPhase.Common.SolutesData.NaOH_dis,
    useLiquidMassFraction=false,
    useLiquidMolality=true,
    refml={1,1,0},
    useGasMassFraction=true,
    refXg={1},
    usePhaseMassFraction=true,
    usePhaseMoleFraction=false));
end NaOH_dis;
