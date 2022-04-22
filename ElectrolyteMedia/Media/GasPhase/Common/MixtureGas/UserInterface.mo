within ElectrolyteMedia.Media.GasPhase.Common.MixtureGas;
record UserInterface "User interface to create a medium"
extends Modelica.Icons.Record;

//Gas phase
  parameter Integer nG = 1 "Number of gas species" annotation(Dialog(tab="Gas phase"));

  Media.Common.Types.GasModel GasModel "Gas equation of state" annotation(choicesAllMatching=true, Dialog(tab="Gas phase"));

  String[nG] gasSubstanceNames = {"gas_1"} "Gas species names" annotation(Dialog(tab="Gas phase"));

  GasDataRecord[nG] datag
      "Select single species data record according to present species" annotation (choicesAllMatching=true, Dialog(
        tab="Gas phase",
        group="Single species gas data"));

  InteractionDataRecord interactiong
      "Interaction data record according to present species" annotation (choicesAllMatching=true, Dialog(
        enable=GasModel == Media.Common.Types.GasModel.PengRobinson,
        tab="Gas phase",
        group="Mixture gas data"));

  Modelica.Media.Interfaces.Types.IdealGas.FluidConstants[nG] fluidconstants "Select species fluid property data record according to present species"
                                                                                                                                                     annotation (choicesAllMatching=true, Dialog(
        tab="Gas phase",
        group="Flui species gas data"));

//General
  String mediumName "Name of medium written as a string" annotation(Dialog(tab="General"));

  SI.Temperature Tstart "Start temperature of medium" annotation(Dialog(tab="General"));

  SI.Pressure pstart "Start pressure of medium" annotation(Dialog(tab="General"));

  Boolean useMassFraction = false "Check box if mass fraction is used" annotation(choices(checkBox = true),Dialog(enable = (not useMoleFraction),tab="General",group = "Reference composition"));

  Boolean useMoleFraction = false "Check box if mole fraction is used" annotation(choices(checkBox = true),Dialog(enable = (not useMassFraction),tab="General",group = "Reference composition"));

  SI.MassFraction[nG] refXg = fill(1/nG,nG) "Reference mass fraction of gas species if useMassFraction" annotation(Dialog(enable = useMassFraction, tab="General",group = "Reference composition"));

  SI.MassFraction[nG] refYg = fill(1/nG,nG) "Reference mole fraction of gas species if useMoleFraction" annotation(Dialog(enable = useMoleFraction, tab="General",group = "Reference composition"));

  final SI.MassFraction[nG] refX = if useMassFraction then refXg elseif useMoleFraction then moleToMassFractions(refYg,MMX[1:nG]) else fill(1/nG,nG);

end UserInterface;
