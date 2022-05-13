within ElectrolyteMedia.Media.SolidLiquidPhase.Common.MixtureLiquid;
record UserInterface
  extends Modelica.Icons.Record;

  //Solid phase
  parameter Integer ns = 1 "Number of solid species" annotation(Dialog(tab="Solid phase"));

  String[:] solidSubstanceNames = {"solid_1"} "Solid species names" annotation(Dialog(tab="Solid phase"));

  DataRecordS[ns] datas "Species data record" annotation(choicesAllMatching=true, Dialog(tab="Solid phase"));

  //Liquid phase
  parameter Integer nL = 2 "Number of liquid species (solutes+solvent)" annotation(Evaluate = true, Dialog(tab="Liquid phase"));

  Media.Common.Types.LiquidModel LiquidModel "Liquid Gibbs excess model" annotation(choicesAllMatching=true, Dialog(tab="Liquid phase"));

  String[:] liquidSubstanceNames = {"solute","solvent"} "Liquid species names (solutes+solvent)" annotation(Evaluate = true,Dialog(tab="Liquid phase"));

  DataRecordL[nL-1] datal "Solutes data record" annotation(choicesAllMatching=true, Dialog(tab="Liquid phase"));

  LiquidInteractionDataRecord interactionL "Binary interaction data record" annotation (choicesAllMatching=true, Dialog(tab="Liquid phase"));

  //Reaction
  parameter Integer nR = 1 "Number of solid-liquid and dissociation equilibria" annotation(Dialog(tab="Reaction"));

  parameter Real[nR,ns+nL] nu = cat(2,identity(nR),-ones(nR,ns+nL-nR)) "Stoichiometry matrix of solid-liquid and dissociation equilibria" annotation(Dialog(tab="Reaction"));

  //General
  String mediumName = "mediumName" "Name of medium" annotation(Dialog(tab="General"));

  SI.Temperature Tstart = 298.15 "Start temperature of medium" annotation(Dialog(tab="General"));

  SI.Pressure pstart = 1e5 "Start pressure of medium" annotation(Dialog(tab="General"));

  Boolean useLiquidMassFraction = false "check if mass fraction for liquid species" annotation(choices(checkBox = true),Dialog(enable = (not useLiquidMoleFraction and not useLiquidMolality),tab="General",group = "Reference composition"));

  Boolean useLiquidMoleFraction = false "check if molality for liquid species" annotation(choices(checkBox = true),Dialog(enable = (not useLiquidMassFraction and not useLiquidMolality),tab="General",group = "Reference composition"));

  Boolean useLiquidMolality = false "check if molality for liquid species" annotation(choices(checkBox = true),Dialog(enable = (not useLiquidMoleFraction and not useLiquidMassFraction),tab="General",group = "Reference composition"));

  SI.MassFraction[nL] refXl = fill(1/nL,nL) "Reference mass fraction of liquid species (solutes+solvent) if useLiquidMassFraction" annotation(Dialog(enable = useLiquidMassFraction, tab="General",group = "Reference composition"));

  SI.MassFraction[nL] refYl = fill(1/nL,nL) "Reference mole fraction of liquid species (solutes+solvent) if useLiquidMoleFraction" annotation(Dialog(enable = useLiquidMoleFraction, tab="General",group = "Reference composition"));

  SI.MassFraction[nL-1] refml = fill(1e-5,nL-1) "Reference molality of liquid species (solutes) if useLiquidMolality" annotation(Dialog(enable = useLiquidMolality, tab="General",group = "Reference composition"));

  Boolean useSolidMassFraction = false "check if mass fraction for solid species" annotation(choices(checkBox = true),Dialog(enable = (not useSolidMoleFraction and not useSolidVolumeFraction),tab="General",group = "Reference composition"));

  Boolean useSolidMoleFraction = false "check if molality for solid species" annotation(choices(checkBox = true),Dialog(enable = (not useSolidMassFraction and not useSolidVolumeFraction),tab="General",group = "Reference composition"));

  Boolean useSolidVolumeFraction = false "check if mass fraction for solid species" annotation(choices(checkBox = true),Dialog(enable = (not useSolidMoleFraction and not useSolidMassFraction),tab="General",group = "Reference composition"));

  SI.MassFraction[ns] refXs = fill(1/ns,ns) "Reference mass fraction of solid species if useSolidMassFraction" annotation(Dialog(enable = useSolidMassFraction, tab="General",group = "Reference composition"));

  SI.MassFraction[ns] refYs = fill(1/ns,ns) "Reference mole fraction of solid species if useSolidMoleFraction" annotation(Dialog(enable = useSolidMoleFraction, tab="General",group = "Reference composition"));

  SI.MassFraction[ns] refPhis = fill(1/ns,ns) "Reference volume fraction of solid species if useSolidMoleFraction" annotation(Dialog(enable = useSolidVolumeFraction,  tab="General",group = "Reference composition"));

  Boolean usePhaseMassFraction = false "Check if liquid phase mass fraction" annotation(choices(checkBox = true),Dialog(enable = (not usePhaseMoleFraction and not usePhaseVolumeFraction),tab="General",group = "Reference composition"));

  Boolean usePhaseMoleFraction = false "Check if liquid phase mole fraction" annotation(choices(checkBox = true),Dialog(enable = (not usePhaseMassFraction and not usePhaseVolumeFraction),tab="General",group = "Reference composition"));

  Boolean usePhaseVolumeFraction = false "Check if liquid phase volume fraction" annotation(choices(checkBox = true),Dialog(enable = (not usePhaseMoleFraction and not usePhaseMassFraction),tab="General",group = "Reference composition"));

  Real refLiquidMassFraction = 0.5 "Reference liquid phase mass fraction" annotation(Dialog(enable = usePhaseMassFraction,tab="General",group = "Reference composition"));

  Real refLiquidMoleFraction = 0.5 "Reference liquid phase mole fraction" annotation(Dialog(enable = usePhaseMoleFraction,tab="General",group = "Reference composition"));

  Real refLiquidVolumeFraction = 0.5 "Reference liquid phase mass fraction" annotation(Dialog(enable = usePhaseVolumeFraction,tab="General",group = "Reference composition"));

  final SI.MassFraction[nL] Xl = if useLiquidMassFraction then refXl elseif useLiquidMoleFraction then moleToMassFractions(refYl,MMX[1+ns:ns+nL]) elseif useLiquidMolality then Functions.LiquidFunctions.calc_Xfromm(refml) else fill(1/nL,nL);

  final SI.MassFraction[ns] Xs = if useSolidMassFraction then refXs elseif useSolidMoleFraction then moleToMassFractions(refYs,MMX[1:ns]) elseif useSolidVolumeFraction then Functions.SolidFunctions.calc_XfromPhi(Tstart,pstart,refPhis) else fill(1/ns,ns);

  final SI.SpecificVolume vl = Functions.LiquidFunctions.calc_v(Tstart,pstart,Xl);

  final SI.SpecificVolume vs = Functions.SolidFunctions.calc_v(Tstart,pstart,Xs);

  final Real LiquidMassFractionFromVolume = refLiquidVolumeFraction*vs/((1-refLiquidVolumeFraction)*vl + refLiquidVolumeFraction*vs);

  final SI.MoleFraction[nL] Yl = massToMoleFractions(Xl,MMX[1+ns:ns+nL]);

  final SI.MoleFraction[ns] Ys = massToMoleFractions(Xs,MMX[1:ns]);

  final SI.MassFraction[ns+nL] refXfull = if usePhaseMassFraction then cat(1,(1-refLiquidMassFraction)*Xs,refLiquidMassFraction*Xl) elseif usePhaseMoleFraction then moleToMassFractions(cat(1,(1-refLiquidMoleFraction)*Ys,refLiquidMoleFraction*Yl),MMX)                                                                                                                                                                                                         elseif usePhaseVolumeFraction then cat(1,(1-LiquidMassFractionFromVolume)*Xs,LiquidMassFractionFromVolume*Xl) else fill(1/(ns+nL),ns+nL);

  final SI.MassFraction[nX] refX= calc_Xred(refXfull);

end UserInterface;
