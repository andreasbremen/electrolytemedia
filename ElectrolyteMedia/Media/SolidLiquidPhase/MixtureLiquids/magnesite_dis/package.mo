within ElectrolyteMedia.Media.SolidLiquidPhase.MixtureLiquids;
package magnesite_dis
extends Common.MixtureLiquid(userInterface(
      ns=1,
      solidSubstanceNames={"MgCO3"},
      datas=SolidPhase.Common.SolidData.magnesite,
      nL=5,
    LiquidModel=ElectrolyteMedia.Media.Common.Types.LiquidModel.DebyeHueckel,
      liquidSubstanceNames={"Mg+2","CO3-2","Na+","Cl-","H2O"},
      datal=LiquidPhase.Common.SolutesData.magnesite_H2O,
      interactionL=LiquidPhase.Common.MixtureSolutesData.magnesite_H2O,
      nR=1,
      nu=[-1,1,1,0,0,0],
      mediumName="Magnesite",
      Tstart=298.15,
      pstart=100000,
    useLiquidMassFraction=false,
    useLiquidMoleFraction=false,
    useLiquidMolality=true,
      refml={0,0,1,1},
      useSolidMassFraction=true,
      usePhaseMassFraction=true,
    refLiquidMassFraction=0.5));
end magnesite_dis;
