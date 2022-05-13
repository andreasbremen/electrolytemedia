within ElectrolyteMedia.Media.SolidLiquidPhase.MixtureLiquids;
package NaCl
extends Common.MixtureLiquid(userInterface(
      ns=1,
      solidSubstanceNames={"NaCl"},
      datas=SolidPhase.Common.SolidData.halite,
    nL=3,
    LiquidModel=ElectrolyteMedia.Media.Common.Types.LiquidModel.DebyeHueckel,
      liquidSubstanceNames={"Na+","Cl-","H2O"},
      datal=LiquidPhase.Common.SolutesData.NaCl,
    interactionL=LiquidPhase.Common.MixtureSolutesData.NaCl,
    nR=1,
    nu=[-1,1,1,0],
      mediumName="NaCl",
    Tstart=298.15,
    pstart=100000,
    useLiquidMassFraction=false,
    useLiquidMoleFraction=false,
    useLiquidMolality=true,
    refml={1,1},
      useSolidMassFraction=true,
      usePhaseMassFraction=true,
    refLiquidMassFraction=0.5));
end NaCl;
