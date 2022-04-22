within ElectrolyteMedia.Media.GasLiquidPhase.Common.MixtureLiquid;
record UserInterface
  extends Modelica.Icons.Record;

  //Gas phase
  parameter Integer nG = 1 "Number of gas species" annotation(Dialog(tab="Gas phase"));

  Media.Common.Types.GasModel GasModel "Gas equation of state" annotation(choicesAllMatching=true, Dialog(tab="Gas phase"));

  String[nG] gasSubstanceNames = {"gas_1"} "Gas species names" annotation(Dialog(tab="Gas phase"));

  Media.GasPhase.Common.GasDataRecord[nG] datag "Select single species data record according to present species" annotation (choicesAllMatching=true, Dialog(
      tab="Gas phase",
      group="Single species gas data"));

  Media.GasPhase.Common.InteractionDataRecord interactiong "Interaction data record according to present species" annotation (choicesAllMatching=true, Dialog(
      enable=GasModel == Media.Common.Types.GasModel.PengRobinson,
      tab="Gas phase",
      group="Mixture gas data"));

  //Liquid phase
  parameter Integer nL = 2 "Number of liquid species (solutes+solvent)" annotation(Evaluate = true, Dialog(tab="Liquid phase"));

  Media.Common.Types.LiquidModel LiquidModel "Liquid Gibbs excess model" annotation(choicesAllMatching=true, Dialog(tab="Liquid phase"));

  String[nL] liquidSubstanceNames = {"solute","solvent"} "Liquid species names (solutes+solvent)" annotation(Evaluate = true,Dialog(tab="Liquid phase"));

  Media.LiquidPhase.Common.DataRecord[nL-1] datal "Solutes data record" annotation(choicesAllMatching=true, Dialog(tab="Liquid phase"));

  Media.LiquidPhase.Common.LiquidInteractionDataRecord interactionl(nLi = nL-1) "Binary interaction data record" annotation (choicesAllMatching=true, Dialog(tab="Liquid phase"));

  //Reaction
  parameter Integer nR_GL = 1 "Number of gas-liquid equilibria" annotation(Dialog(tab="Reaction"));

  Real[nR_GL,nG+nL] nu_GL = zeros(nR_GL,nG+nL) "Stoichiometry matrix of gas-liquid and dissociation equilibria" annotation(Dialog(tab="Reaction"));

  parameter Integer nR_L = 1 "Number of dissociation equilibria" annotation(Dialog(tab="Reaction"));

  Real[nR_L,nL] nu_L = zeros(nR_L,nL) "Stoichiometry matrix of dissociation equilibria" annotation(Dialog(tab="Reaction"));

  final parameter Integer nR = nR_GL + nR_L "Number of gas-liquid and dissociation equilibria";

  final Real[nR,nG+nL] nu = {if i<nR_GL+1 then nu_GL[i,:] else cat(1,zeros(nG),nu_L[i-nR_GL,:]) for i in 1:nR} "Stoichiometry matrix of gas-liquid and dissociation equilibria" annotation(Dialog(tab="Reaction"));

  //General
  String mediumName "Name of medium written as a string" annotation(Dialog(tab="General"));

  SI.Temperature Tstart "Start temperature of medium" annotation(Dialog(tab="General"));

  SI.Pressure pstart "Start pressure of medium" annotation(Dialog(tab="General"));

  Boolean useLiquidMassFraction = false "Check box if mass fraction is used for liquid species" annotation(choices(checkBox = true),Dialog(enable = (not useLiquidMoleFraction and not useLiquidMolality),tab="General",group = "Reference composition"));

  Boolean useLiquidMoleFraction = false "Check box if mole fraction is used for liquid species" annotation(choices(checkBox = true),Dialog(enable = (not useLiquidMassFraction and not useLiquidMolality),tab="General",group = "Reference composition"));

  Boolean useLiquidMolality = false "Check box if molality is used for liquid species" annotation(choices(checkBox = true),Dialog(enable = (not useLiquidMoleFraction and not useLiquidMassFraction),tab="General",group = "Reference composition"));

  SI.MassFraction[nL] refXl = fill(1/nL,nL) "Reference mass fraction of liquid species (solutes+solvent) if useLiquidMassFraction" annotation(Dialog(enable = useLiquidMassFraction, tab="General",group = "Reference composition"));

  SI.MassFraction[nL] refYl = fill(1/nL,nL) "Reference mole fraction of liquid species (solutes+solvent) if useLiquidMoleFraction" annotation(Dialog(enable = useLiquidMoleFraction, tab="General",group = "Reference composition"));

  SI.MassFraction[nL-1] refml = fill(1e-5,nL-1) "Reference molality of liquid species (solutes) if useLiquidMolality" annotation(Dialog(enable = useLiquidMolality, tab="General",group = "Reference composition"));

  Boolean useGasMassFraction = false "Check box if mass fraction is used for gaseous species" annotation(choices(checkBox = true),Dialog(enable = (not useGasMoleFraction),tab="General",group = "Reference composition"));

  Boolean useGasMoleFraction = false "Check box if mole fraction is used for gaseous species" annotation(choices(checkBox = true),Dialog(enable = (not useGasMassFraction),tab="General",group = "Reference composition"));

  SI.MassFraction[nG] refXg = fill(1/nG,nG) "Reference mass fraction of gas species if useGasMassFraction" annotation(Dialog(enable = useGasMassFraction, tab="General",group = "Reference composition"));

  SI.MassFraction[nG] refYg = fill(1/nG,nG) "Reference mole fraction of gas species if useGasMoleFraction" annotation(Dialog(enable = useGasMoleFraction, tab="General",group = "Reference composition"));

  Boolean usePhaseMassFraction = false "Check box if mass fraction is used for liquid phase" annotation(choices(checkBox = true),Dialog(enable = (not usePhaseMoleFraction),tab="General",group = "Reference composition"));

  Boolean usePhaseMoleFraction = false "Check box if mole fraction is used for liquid phase" annotation(choices(checkBox = true),Dialog(enable = (not usePhaseMassFraction),tab="General",group = "Reference composition"));

  Real refLiquidMassFraction = 0.5 "Reference liquid phase mass fraction" annotation(Dialog(enable = usePhaseMassFraction,tab="General",group = "Reference composition"));

  Real refLiquidMoleFraction = 0.5 "Reference liquid phase mole fraction" annotation(Dialog(enable = usePhaseMoleFraction,tab="General",group = "Reference composition"));

  final SI.MassFraction[nL] Xl = if useLiquidMassFraction then refXl elseif useLiquidMoleFraction then moleToMassFractions(refYl,MMX[1+nG:nG+nL]) elseif useLiquidMolality then Functions.LiquidFunctions.calc_Xfromm(refml) else fill(1/nL,nL);

  final SI.MassFraction[nG] Xg = if useGasMassFraction then refXg elseif useGasMoleFraction then moleToMassFractions(refYg,MMX[1:nG]) else fill(1/nG,nG);

  final SI.MoleFraction[nL] Yl = massToMoleFractions(Xl,MMX[1+nG:nG+nL]);

  final SI.MoleFraction[nG] Yg = massToMoleFractions(Xg,MMX[1:nG]);

  final SI.MassFraction[nG+nL] refXfull = if usePhaseMassFraction then cat(1,(1-refLiquidMassFraction)*Xg,refLiquidMassFraction*Xl) elseif usePhaseMoleFraction then moleToMassFractions(cat(1,(1-refLiquidMoleFraction)*Yg,refLiquidMoleFraction*Yl),MMX) else fill(1/(nG+nL),nG+nL);

  final SI.MassFraction[nX] refX = calc_Xred(refXfull);

end UserInterface;
