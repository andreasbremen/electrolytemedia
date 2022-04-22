within ElectrolyteMedia.Media.LiquidPhase.MixtureLiquids;
package AceticAcid_NaOH
extends Common.MixtureLiquid(userInterface(
    nL=7,
    LiquidModel=ElectrolyteMedia.Media.Common.Types.LiquidModel.Ideal,
    substanceNames={"AceticAcid","Acetate","Sodium","Chloride","Hydroxide",
        "HydrogenIon","Water"},
    data=Common.SolutesData.Acetic_NaCl,
    interaction=Common.MixtureSolutesData.Acetic_NaCl,
    nR=2,
    nu=[-1,1,0,0,0,1,0; 0,0,0,0,1,1,-1],
    mediumName="AceticAcidWithNaCl",
    Tstart=298.15,
    pstart=100000,
    useLiquidMolality=true,
    refml={1,0,1,1,0,0}));                                                //Ideal

end AceticAcid_NaOH;
