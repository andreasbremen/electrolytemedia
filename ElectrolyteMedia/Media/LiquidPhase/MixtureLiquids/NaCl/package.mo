within ElectrolyteMedia.Media.LiquidPhase.MixtureLiquids;
package NaCl
extends Common.MixtureLiquid(userInterface(
    nL=3,
    LiquidModel=ElectrolyteMedia.Media.Common.Types.LiquidModel.Pitzer,
    substanceNames={"Sodium","Chloride","Water"},
    data=Common.SolutesData.NaCl,
    interaction=Common.MixtureSolutesData.NaCl,
    nR=0,
    mediumName="NaCl",
    Tstart=298.15,
    pstart=100000,
    useLiquidMolality=true,
    refml={1,1}));

end NaCl;
