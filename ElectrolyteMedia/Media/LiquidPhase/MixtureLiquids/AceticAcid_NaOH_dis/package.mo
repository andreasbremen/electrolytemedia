within ElectrolyteMedia.Media.LiquidPhase.MixtureLiquids;
package AceticAcid_NaOH_dis
extends Common.MixtureLiquid(userInterface(
      nL=8,
    LiquidModel=ElectrolyteMedia.Media.Common.Types.LiquidModel.Ideal,
    substanceNames={"AceticAcid","Acetate","SodiumHydroxide","Sodium","Chloride","Hydroxide",
        "HydrogenIon","Water"},
      data=Common.SolutesData.Acetic_NaOH,
      interaction=Common.MixtureSolutesData.Acetic_NaOH,
      nR=3,
      nu=[-1,1,0,0,0,0,1,0; 0,0,0,0,0,1,1,-1; 0,0,-1,1,0,1,0,0],
    mediumName="AceticAcidWithNaCl",
    Tstart=298.15,
    pstart=100000,
    useLiquidMolality=true,
    refml={1,0,1,0,0,0,0}));

end AceticAcid_NaOH_dis;
