within ElectrolyteMedia.Media.SolidLiquidPhase.MixtureLiquids;
package ExampleMedium
extends Common.MixtureLiquid(userInterface(
      ns=1,
      solidSubstanceNames={"NaCl"},
      datas=SolidPhase.Common.SolidData.halite,
      nL=5,
    LiquidModel=ElectrolyteMedia.Media.Common.Types.LiquidModel.Bromley,
      liquidSubstanceNames={"Na+","Cl-","OH-","H+","H2O"},
      datal=LiquidPhase.Common.SolutesData.NaCl_H2O,
    interactionL=LiquidPhase.Common.MixtureSolutesData.NaCl_H2O,
    nR=1,
    nu=[-1,1,1,0,0,0],
      mediumName="NaCl",
    Tstart=298.15,
    pstart=100000,
    useLiquidMassFraction=false,
    useLiquidMoleFraction=true,
    useLiquidMolality=false,
    refYl={0,0,1e-7,1e-7,1 - 2e-7},
      useSolidMassFraction=true,
      usePhaseMassFraction=true,
    refLiquidMassFraction=0.5));
end ExampleMedium;
