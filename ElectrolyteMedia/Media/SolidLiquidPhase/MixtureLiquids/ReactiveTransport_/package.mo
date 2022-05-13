within ElectrolyteMedia.Media.SolidLiquidPhase.MixtureLiquids;
package ReactiveTransport_
extends Common.MixtureLiquid(userInterface(
      ns=3,
      solidSubstanceNames={"Quartz","Calcite","Dolomite"},
      datas=SolidPhase.Common.SolidData.ReactiveTransport,
      nL=11,
    LiquidModel=ElectrolyteMedia.Media.Common.Types.LiquidModel.Ideal,
      liquidSubstanceNames={"CO2","OH-","HCO3-","CO3-2","Mg+2","Ca+2","Na+","Cl-","SiO2","H+",
          "H2O"},
      datal=LiquidPhase.Common.SolutesData.ReactiveTransport_,
      interactionL=LiquidPhase.Common.MixtureSolutesData.ReactiveTransport_,
      nR=6,
      nu=[1,0,0,0,0,0,0,0,0,0,0,-1,0,0; 0,1,0,0,0,0,-1,0,-1,0,0,0,0,0; 0,0,1,0,0,
          0,-2,-1,-1,0,0,0,0,0; 0,0,0,1,0,0,-1,0,0,0,0,0,-2,1; 0,0,0,0,1,0,-1,0,
          0,0,0,0,-1,0; 0,0,0,0,0,1,0,0,0,0,0,0,1,-1],
      mediumName="ReactiveTransport",
      Tstart=333.15,
      pstart=1000000,
    useLiquidMassFraction=false,
    useLiquidMoleFraction=false,
    useLiquidMolality=true,
    refml={0,0,0,0,0,0,0,0,0,0},
    useSolidMassFraction=false,
    useSolidVolumeFraction=true,
    refPhis={0.98,0.02-1e-5,1e-5},
    usePhaseMassFraction=false,
    usePhaseVolumeFraction=true,
    refLiquidMassFraction=1,
    refLiquidVolumeFraction=0.1));
                                //{0.75,0,0,0.05,0.01,0.9,1.02,0,0,0},
end ReactiveTransport_;