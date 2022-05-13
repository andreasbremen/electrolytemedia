within ElectrolyteMedia.Media.SolidLiquidPhase.MixtureLiquids;
package NaCl_dis
extends Common.MixtureLiquid(userInterface(
      ns=1,
      solidSubstanceNames={"NaCl"},
      datas=SolidPhase.Common.SolidData.halite,
      nL=5,
    LiquidModel=ElectrolyteMedia.Media.Common.Types.LiquidModel.Ideal,
      liquidSubstanceNames={"Na+","Cl-","OH-","H+","H2O"},
      datal=LiquidPhase.Common.SolutesData.NaCl_H2O,
      interactionL=LiquidPhase.Common.MixtureSolutesData.NaCl_H2O,
      nR=2,
      nu=[-1,1,1,0,0,0;0,0,0,1,1,-1],
      mediumName="NaCl",
    Tstart=298.15,
    pstart=100000,
    useLiquidMassFraction=false,
    useLiquidMoleFraction=false,
    useLiquidMolality=true,
      refml={1,1,0,0},
      useSolidMassFraction=true,
      usePhaseMassFraction=true,
    refLiquidMassFraction=0.5));
end NaCl_dis;
