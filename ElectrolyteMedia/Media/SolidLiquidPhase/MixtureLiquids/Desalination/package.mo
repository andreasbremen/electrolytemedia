within ElectrolyteMedia.Media.SolidLiquidPhase.MixtureLiquids;
package Desalination
extends Common.MixtureLiquid(userInterface(
      ns=3,
      solidSubstanceNames={"NaCl","CaSO4.2H2O","SiO2"},
      datas=SolidPhase.Common.SolidData.Desal_1,
      nL=10,
      LiquidModel=ElectrolyteMedia.Media.Common.Types.LiquidModel.Ideal,
      liquidSubstanceNames={"Na+","K+","Ca+2","Mg+2","Cl-","SO4-2","SiO2","OH-","H+","H2O"},
      datal=LiquidPhase.Common.SolutesData.Seawater,
      interactionL=LiquidPhase.Common.MixtureSolutesData.Seawater,
      nR=4,
      nu=[0,0,0,0,0,0,0,0,0,0,1,1,-1; -1,0,0,1,0,0,0,1,0,0,0,0,0; 0,-1,0,0,0,1,0,
          0,1,0,0,0,2; 0,0,-1,0,0,0,0,0,0,1,0,0,0],
      mediumName="Desalination",
      Tstart=298.15,
      pstart=100000,
    useLiquidMassFraction=false,
      useLiquidMoleFraction=false,
      useLiquidMolality=true,
      refml={0.1,0,0.1,0,0.1,0.1,0.01,0,0},
      useSolidMassFraction=true,
      usePhaseMassFraction=true,
      refLiquidMassFraction=1));

end Desalination;
