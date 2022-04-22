within ElectrolyteMedia.Media.LiquidPhase.MixtureLiquids;
package PhosphoricAcid_NaOH
extends Common.MixtureLiquid(userInterface(
    nL=8,
    LiquidModel=ElectrolyteMedia.Media.Common.Types.LiquidModel.DebyeHueckel,
    substanceNames={"PhosphoricAcid","DihydrogenPhosphate","HydrogenPhosphate",
        "Phosphate","Sodium","Hydroxide","HydrogenIon","Water"},
    data=Common.SolutesData.PhosphoricAcid,
      interaction=Common.MixtureSolutesData.PhosphoricAcid,
    nR=4,
    nu=[-1,1,0,0,0,0,1,0; 0,-1,1,0,0,0,1,0; 0,0,-1,1,0,0,1,0; 0,0,0,0,0,1,1,-1],
    mediumName="PhosphoricAcidwithNaOH",
    Tstart=298.15,
    pstart=100000,
    useLiquidMolality=true,
    refml={1,0,0,0,0.1,0.1,0}));


end PhosphoricAcid_NaOH;
