within ElectrolyteMedia.Media.LiquidPhase.Common.MixtureLiquid;
record UserInterface "User interface to create a medium"
  extends Modelica.Icons.Record;

  //Liquid phase
  parameter Integer nL = 2 "Number of liquid species (solutes+solvent)" annotation(Evaluate = true, Dialog(tab="Liquid phase"));

  Media.Common.Types.LiquidModel LiquidModel "Liquid Gibbs excess model" annotation(choicesAllMatching=true, Dialog(tab="Liquid phase"));

  String[nL] substanceNames = {"solute","solvent"} "Liquid species names (solutes+solvent)" annotation(Evaluate = true,Dialog(tab="Liquid phase"));

  Media.LiquidPhase.Common.DataRecord[nL-1] data "Solutes data record" annotation(choicesAllMatching=true, Dialog(tab="Liquid phase"));

  LiquidInteractionDataRecord interaction(nLi = nL-1) "Binary interaction data record"
    annotation (choicesAllMatching=true, Dialog(tab="Liquid phase"));

  //Reaction
  parameter Integer nR = 0 "Number of dissociation equilibria" annotation(Dialog(tab="Reaction"));

  Real[nR,nL] nu = zeros(nR,nL) "Stoichiometry matrix of gas-liquid and dissociation equilibria" annotation(Dialog(tab="Reaction"));

  //General
  String mediumName = "Name of medium" "Name of medium written as a string" annotation(Dialog(tab="General"));

  SI.Temperature Tstart = 298.15 "Start temperature of medium" annotation(Dialog(tab="General"));

  SI.Pressure pstart = 1e5 "Start pressure of medium" annotation(Dialog(tab="General"));

  Boolean useLiquidMassFraction = false "Check box if mass fraction is used for liquid species" annotation(choices(checkBox = true),Dialog(enable = (not useLiquidMoleFraction and not useLiquidMolality),tab="General",group = "Reference composition"));

  Boolean useLiquidMoleFraction = false "Check box if mole fraction is used for liquid species" annotation(choices(checkBox = true),Dialog(enable = (not useLiquidMassFraction and not useLiquidMolality),tab="General",group = "Reference composition"));

  Boolean useLiquidMolality = false "Check box if molality is used for liquid species" annotation(choices(checkBox = true),Dialog(enable = (not useLiquidMoleFraction and not useLiquidMassFraction),tab="General",group = "Reference composition"));

  SI.MassFraction[nL] refXl = fill(1/nL,nL) "Reference mass fraction of liquid species (solutes+solvent) if useLiquidMassFraction" annotation(Dialog(enable = useLiquidMassFraction, tab="General",group = "Reference composition"));

  SI.MassFraction[nL] refYl = fill(1/nL,nL) "Reference mole fraction of liquid species (solutes+solvent) if useLiquidMoleFraction" annotation(Dialog(enable = useLiquidMoleFraction, tab="General",group = "Reference composition"));

  SI.MassFraction[nL-1] refml = fill(1e-5,nL-1) "Reference molality of liquid species (solutes) if useLiquidMolality" annotation(Dialog(enable = useLiquidMolality, tab="General",group = "Reference composition"));

  final SI.MassFraction[nL] refXfull = if useLiquidMassFraction then refXl elseif useLiquidMoleFraction then moleToMassFractions(refYl,MMX) elseif useLiquidMolality then Functions.calc_Xfromm(refml) else fill(1/nL,nL);

  final SI.MassFraction[nX] refX=calc_Xred_from_Xfull(refXfull);

end UserInterface;
