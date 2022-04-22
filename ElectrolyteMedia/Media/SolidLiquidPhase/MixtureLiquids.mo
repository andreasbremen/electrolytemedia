within ElectrolyteMedia.Media.SolidLiquidPhase;
package MixtureLiquids
  extends Modelica.Icons.VariantsPackage;

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

  package KCl_NaCl
  extends Common.MixtureLiquid(userInterface(
        ns=2,
        solidSubstanceNames={"KCl","NaCl"},
        datas=SolidPhase.Common.SolidData.KCl_NaCl,
        nL=6,
        LiquidModel=ElectrolyteMedia.Media.Common.Types.LiquidModel.Ideal,
        liquidSubstanceNames={"K+","Cl-","Na+","OH-","H+","H2O"},
        datal=LiquidPhase.Common.SolutesData.KCl_NaCl,
        interactionL=LiquidPhase.Common.MixtureSolutesData.KCl_NaCl,
      nR=3,
      nu=[-1,0,1,1,0,0,0,0; 0,-1,0,1,1,0,0,0; 0,0,0,0,0,1,1,-1],
        mediumName="NaCl",
      Tstart=298.15,
      pstart=100000,
      useLiquidMassFraction=false,
      useLiquidMoleFraction=true,
      useLiquidMolality=false,
        refYl={0,0,0,1e-7,1e-7,1 - 2e-7},
        useSolidMassFraction=true,
        usePhaseMassFraction=true,
      refLiquidMassFraction=0.5));       //refml={1,2,1,0,0,0},
      //refYl={0,0,0,1e-7,1e-7,1 - 2e-7},
                                  //true,
                             //false,
  end KCl_NaCl;

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

  package NaCl_dis_
  extends Common.MixtureLiquid(userInterface(
        ns=1,
        solidSubstanceNames={"NaCl"},
        datas=SolidPhase.Common.SolidData.halite,
        nL=5,
      LiquidModel=ElectrolyteMedia.Media.Common.Types.LiquidModel.Pitzer,
        liquidSubstanceNames={"Na+","Cl-","OH-","H+","H2O"},
        datal=LiquidPhase.Common.SolutesData.NaCl_H2O,
        interactionL=LiquidPhase.Common.MixtureSolutesData.NaCl_H2O,
        nR=1,
        nu=[-1,1,1,0,0,0],
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
  end NaCl_dis_;

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

  package ReactiveTransport
  extends Common.MixtureLiquid(userInterface(
        ns=3,
        solidSubstanceNames={"Quartz","Calcite","Dolomite"},
        datas=SolidPhase.Common.SolidData.ReactiveTransport,
        nL=11,
      LiquidModel=ElectrolyteMedia.Media.Common.Types.LiquidModel.Ideal,
        liquidSubstanceNames={"CO2","HCO3-","CO3-2","Mg+2","Ca+2","Na+","Cl-","SiO2","OH-","H+",
            "H2O"},
        datal=LiquidPhase.Common.SolutesData.ReactiveTransport,
        interactionL=LiquidPhase.Common.MixtureSolutesData.ReactiveTransport,
        nR=6,
        nu=[1,0,0,0,0,0,0,0,0,0,-1,0,0,0; 0,0,1,0,0,-2,-1,-1,0,0,0,0,0,0; 0,1,0,0,0,-1,0,-1,0,0,0,0,0,0;0,0,0,-1,1,0,0,0,0,0,0,0,1,-1; 0,0,0,0,-1,1,0,0,0,0,
            0,0,1,0; 0,0,0,0,0,0,0,0,0,0,0,-1,-1,1],
        mediumName="ReactiveTransport",
        Tstart=333.15,
        pstart=1000000,
      useLiquidMassFraction=false,
      useLiquidMoleFraction=false,
      useLiquidMolality=true,
      refml={0,0,0,0,0,0,0,0,0,0},
      useSolidMassFraction=false,
      useSolidVolumeFraction=true,
        refPhis={0.98,0.02 - 0,0},
      usePhaseMassFraction=false,
      usePhaseVolumeFraction=true,
      refLiquidVolumeFraction=0.5));
                                  //{0.75,0,0,0.05,0.01,0.9,1.02,0,0,0},
  end ReactiveTransport;

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
end MixtureLiquids;
