within ElectrolyteMedia.Media.LiquidPhase.MixtureLiquids;
package Electrolysis
extends Common.MixtureLiquid(userInterface(
    nL=14,
    LiquidModel=ElectrolyteMedia.Media.Common.Types.LiquidModel.Ideal,
    substanceNames={"1-H+","2-OH-","3-CO","4-CO2","5-K+","6-HCO3-","7-CO3-2","8-KSO4-","9-SO4-2","10-HSO4-","11-KOH","12-H2","13-O2","14-Water"},
    data=Common.SolutesData.Electrolysis,
    interaction=Common.MixtureSolutesData.Electrolysis,
    nR=6,
      nu=[-1,-1,0,0,0,0,0,0,0,0,0,0,0,1; -1,0,0,1,0,-1,0,0,0,0,0,0,0,1; -1,0,
          0,0,0,1,-1,0,0,0,0,0,0,0; 0,1,0,0,1,0,0,0,0,0,-1,0,0,0; 0,0,0,0,1,0,0,
          -1,1,0,0,0,0,0; -1,0,0,0,0,0,0,0,-1,1,0,0,0,0],
    mediumName="Electrolyte for PEM Electrolysis",
      Tstart=298.15,
      pstart=100000,
    useLiquidMolality=true,
      refml={0,0,0,0.0375,1,0,0,0,0.5,0,0.0,0,0}));

end Electrolysis;
