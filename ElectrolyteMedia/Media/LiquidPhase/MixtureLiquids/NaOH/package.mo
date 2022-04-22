within ElectrolyteMedia.Media.LiquidPhase.MixtureLiquids;
package NaOH
extends Common.MixtureLiquid(userInterface(
      nL=4,
    LiquidModel=ElectrolyteMedia.Media.Common.Types.LiquidModel.Pitzer,
      substanceNames={"Sodium","Hydroxide","NaOH","Water"},
      data=Common.SolutesData.NaOH_dis,
      interaction=Common.MixtureSolutesData.NaOH_dis,
      nR=1,
      nu=[-1,-1,1,0],
      mediumName="NaOH",
    Tstart=298.15,
    pstart=100000,
    useLiquidMolality=true,
      refml={0,0,1}));

end NaOH;
