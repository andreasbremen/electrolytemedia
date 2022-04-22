within ElectrolyteMedia.Media.LiquidPhase.MixtureLiquids;
package HydrochloricAcid
extends Common.MixtureLiquid(userInterface(
    nL=6,
    LiquidModel=ElectrolyteMedia.Media.Common.Types.LiquidModel.Pitzer,
      substanceNames={"HydrochloricAcid","Sodium","Chloride","Hydroxide","HydrogenIon","Water"},
      data=Common.SolutesData.HydrochloricAcid,
    interaction=Common.MixtureSolutesData.HydrochloricAcid,
    nR=2,
    nu=[-1,0,1,0,1,0; 0,0,0,1,1,-1],
    mediumName="HydrochloricAcidWithNaOH",
    Tstart=298.15,
    pstart=100000,
    useLiquidMolality=true,
    refml={0.34,0,0,0,0}));

end HydrochloricAcid;
