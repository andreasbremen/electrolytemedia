within ElectrolyteMedia.Media.SolidPhase.Common.MixtureSolid;
record UserInterface "User interface to create new medium"
  extends Modelica.Icons.Record;

  //Solid phase
  parameter Integer nX = 1 "Number of solid species" annotation(Evaluate = true, Dialog(tab="Solid phase"));

  String[nX] substanceNames = {"solid_1"} "Solid species names" annotation(Evaluate = true,Dialog(tab="Solid phase"));

  DataRecord[nX] data "Solid data record"
    annotation (choicesAllMatching=true, Dialog(tab="Solid phase"));

  //General
  String mediumName "Name of medium written as a string" annotation(Dialog(tab="General"));

  SI.Temperature Tstart "Start temperature of medium" annotation(Dialog(tab="General"));

  SI.Pressure pstart "Start pressure of medium" annotation(Dialog(tab="General"));

  Boolean useSolidMassFraction = false "Check box if mass fraction is used for solid species" annotation(choices(checkBox = true),Dialog(enable = (not useSolidMoleFraction),tab="General",group = "Reference composition"));

  Boolean useSolidMoleFraction = false "Check box if mole fraction is used for solid species" annotation(choices(checkBox = true),Dialog(enable = (not useSolidMassFraction),tab="General",group = "Reference composition"));

  SI.MassFraction[nX] refXs = fill(1/nX,nX) "Reference mass fraction of solid species if useSolidMassFraction" annotation(Dialog(enable = useSolidMassFraction, tab="General",group = "Reference composition"));

  SI.MassFraction[nX] refYs = fill(1/nX,nX) "Reference mole fraction of solid species if useSolidMoleFraction" annotation(Dialog(enable = useSolidMoleFraction, tab="General",group = "Reference composition"));

  final SI.MassFraction[nX] refX = if useSolidMassFraction then refXs elseif useSolidMoleFraction then moleToMassFractions(refYs,MMX) else fill(1/nX,nX);

end UserInterface;
