within ElectrolyteMedia.Media.LiquidPhase.MixtureLiquids;
package OxalicAcid_NaOH
extends Common.MixtureLiquid(userInterface(
    nL=7,
    LiquidModel=ElectrolyteMedia.Media.Common.Types.LiquidModel.DebyeHueckel,
    substanceNames={"OxalicAcid","HydrogenOxalate","Oxalate","Sodium",
        "Hydroxide","HydrogenIon","Water"},
    data=Common.SolutesData.OxalicAcid,
    interaction=Common.MixtureSolutesData.OxalicAcid,
    nR=3,
    nu=[-1,1,0,0,0,1,0; 0,-1,1,0,0,1,0; 0,0,0,0,1,1,-1],
    mediumName="OxalicAcidwithNaOH",
    Tstart=298.15,
    pstart=100000,
    useLiquidMolality=true,
    refml={1,0,0,1,1,0}));

end OxalicAcid_NaOH;
