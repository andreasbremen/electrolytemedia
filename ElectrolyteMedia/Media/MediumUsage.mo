within ElectrolyteMedia.Media;
class MediumUsage "Medium Usage"
  extends Modelica.Icons.Information;
  annotation (Documentation(info="<html>
<h4>General information</h4>
<p><br>The <i>Media</i> package is used to define a medium for the use in custom unit models or models from the <i><a href=\"Modelica.Fluid\">Fluid</a></i>  library. A medium allows the user to specify the chemical species, a stoichiometry matrix to account for equilibrium reactions, thermodynamic models, and reference composition to use for a simulation. </p>
<p><br><br>The package <i>Media</i> is based on the modelica library <a href=\"Modelica.Media\">Modelica.Media</a> (see here for basic information). We provide thermodynamic models in combination with isopotential conditions and law of mass action (LMA) to describe phase equilibrium and dissociation equilibrium for phases relevant in aqueous electrolyte thermodynamics, as given in the following:</p>
<table cellspacing=\"0\" cellpadding=\"2\" border=\"1\" width=\"100%\"><tr>
<td><h4>Phase</h4></td>
<td><h4>Model</h4></td>
<td><h4>Path</h4></td>
</tr>
<tr>
<td><p>Gas</p></td>
<td><p>Peng-Robinson equation of state</p><p>Ideal gas law</p></td>
<td><p><a href=\"Media.GasPhase\">Media.GasPhase</a></p></td>
</tr>
<tr>
<td><p>Solid</p></td>
<td><p>Solid phase thermodynamics based on Holland and Powell 2011</p></td>
<td><p><a href=\"Media.SolidPhase\">Media.SolidPhase</a></p></td>
</tr>
<tr>
<td><p>Liquid</p></td>
<td><p>Ideal liquid phase at infinite dilution</p><p>Debye limiting law for single charged ion pairs</p><p>Debye-Hueckel model</p><p>Bromley model</p><p>Pitzer model</p></td>
<td><p><a href=\"Media.LiquidPhase\">Media.LiquidPhase</a></p></td>
</tr>
<tr>
<td><p>Gas and liquid</p></td>
<td><p>Combination of gas and liquid phase models for two phase Media</p></td>
<td><p><a href=\"Media.GasLiquidPhase\">Media.GasLiquidPhase</a></p></td>
</tr>
<tr>
<td><p>Solid and liquid with a variable number of solid phases</p></td>
<td><p>Combination of solid and liquid phase models for multiple phase Media</p></td>
<td><p><a href=\"Media.SolidLiquidPhase\">Media.SolidLiquidPhase</a></p></td>
</tr>
</table>
<p><br>Table 1: Implemented media packages for aqueous electrolyte systems </p>
<h4>Medium incorporating dissocation and phase equilibria</h4>
<p>The main feature of the <i>ElectrolyteMedia</i> library is the simultaneous incorporation of the law of mass action and isopotential conditions in the <i>BaseProperties</i> model to represent dissociation and phase equilibrium, respectively. The number of considered dissociation and phase equilibria reduces the number of independent mass fractions, as they become dependent on each other. As a result, the number of degrees of freedom is reduced, and hence, the independent mass fractions may not be set arbitrarily. Instead, we project the mass fraction vector of the <i>BaseProperties</i> model to the reaction invariants that remain constant at a predefined overall composition. The reaction invariants and two furhter intensive thermodynamic properties are then the degrees of freedom of the <i>BaseProperties</i> model. Mathematically, reaction invariants may be interpreted as a linear combination of the mass-based atom balance.</p>
<p>To calulate the reaction invariants, a mass-based atom balance could be used. In that case, only the atomic composition of considered species is required. With the mass of each kind of atom that compose each species, an atom balance can be written for each kind of atom. However, this renders all species to be in equilibrium. There are, however, applications, where we consider only partial equilibrium of speceies present in a system. Hence, we use a stoichiometric approach, where a mass-based stoichiometry matrix <b><img src=\"modelica://ElectrolyteMedia/Resources/Images/equations/equation-n3oAAX7n.png\" alt=\"nu\"/></b> is defined for a given number of species. The stoichiometry matrix is of size <i>nF </i>x <i>nR</i>, where <i>nF</i> is the full vector of present species and <i>nR</i> is the number of considered equilibrium reactions.</p>
<p>The reaction invariants are then calculated from the nullspace of the stoichiometry matrix <img src=\"modelica://ElectrolyteMedia/Resources/Images/equations/equation-5frUFgXO.png\" alt=\"nu\"/> for which holds:</p>
<p style=\"margin-left: 30px;\"><img src=\"modelica://ElectrolyteMedia/Resources/Images/equations/equation-49jw5Xzy.png\" alt=\"nu * lambda = 0\"/></p>
<p>The nullspace <img src=\"modelica://ElectrolyteMedia/Resources/Images/equations/equation-AO7bzKUX.png\" alt=\"lambda\"/> is calculated from a Gauss-Jordan algorithm provided in <a href=\"modelica://ElectrolyteMedia/Media/Common/Reaction/calc_lambda.mo\">calc_lambda</a>. This function ensures that element of the reaction invariant is greater than zero for the purpose of normalizing <img src=\"modelica://ElectrolyteMedia/Resources/Images/equations/equation-giWA1RiN.png\" alt=\"X\"/> and <img src=\"modelica://ElectrolyteMedia/Resources/Images/equations/equation-TxKB1tUY.png\" alt=\"X_full
\"/>.This allows the use of the reaction variant vector <i>X</i> to be used as mass fraction vector within a <i>BaseProperties</i> model, i.e., they are equivalent. For that purpose, we calculate the reaction invariant from the full vector of mass fractions <img src=\"modelica://ElectrolyteMedia/Resources/Images/equations/equation-XtYIWu7B.png\" alt=\"X_full\"/> by: </p>
<p style=\"margin-left: 30px;\"><img src=\"modelica://ElectrolyteMedia/Resources/Images/equations/equation-VlfPgBJU.png\" alt=\"X = lambda^T*X_full/(sum(lambda^T*X_full))\"/></p>
<p>To distinguish the mass fraction vector X from the independent mass fraction vector <img src=\"modelica://ElectrolyteMedia/Resources/Images/equations/equation-pUQ9bZ6A.png\" alt=\"X_full\"/>, we use the name convention <img src=\"modelica://ElectrolyteMedia/Resources/Images/equations/equation-bjimK7ur.png\" alt=\"X_red\"/> as the mass fraction vector denoting the reaction invariant, i.e., <img src=\"modelica://ElectrolyteMedia/Resources/Images/equations/equation-I9ChVnme.png\" alt=\"X_red = X\"/> with subscirpt <i>red</i> denoting &quot;reduced&quot;.</p>
<h4>Creating a medium</h4>
<p><br>An example of a medium implementation can be found <a href=\"Tests.Models_mabo.Titration.partial_Titration\">here</a>.</p>
<p>To create and use these new media a user interface is provided. This guide only includes information on the user interface. For more information, especially the underlying code, see Modelica&apos;s documentation on <a href=\"Modelica.Media.UsersGuide.MediumDefinition\">Medium definition</a>. </p>
<p><br><img src=\"modelica://ElectrolyteMedia/Resources/Images/Pictures/Extends.PNG\"/></p>
<p>Figure 1: A Medium </p>
<p><br><br>To create a medium, a new package has to be created as seen in Figure 1. This package needs to extend from one of the paths shown in Table 1. After creating the package the diagram view must be enabled to make the user interface visible. The user interface can be opened by double clicking, all necessary information should be filled in by the user. The information to be filled in is shortly explained in comments in the UI. This documentation lists further details if needed. </p>
<p>Since the user interfaces for different phases differ only slightly, only the liquid phase user interface is discussed here.</p>
<p><br><b>Tab: General</b></p>
<p style=\"margin-left: 30px;\">Here, parameters for initialization of the medium are set. In detail, the medium requires a temperature, a pressure, and a composition to calculate initial values for the equilibrium composition.</p>
<p style=\"margin-left: 30px;\"><b>Parameters</b>:</p>
<p style=\"margin-left: 60px;\"><b>Tstart</b>: start temperature. The medium uses this temperature to initialize the reaction equilibria of considered reactions.</p>
<p style=\"margin-left: 60px;\"><b>pstart</b>: start pressure. The medium uses this temperature to initialize the reaction equilibria of considered reactions.</p>
<p style=\"margin-left: 30px;\"><b>Reference composition</b>:</p>
<p style=\"margin-left: 60px;\">The medium uses the reference composition to initialize the reaction equilibria of considered reactions. The reference composition is composed of the undissociated species that are initially present in the considered medium</p>
<p style=\"margin-left: 60px;\"><b>Check boxes</b>: For multiple substance systems. Check box to specify in which way the composition is given. Either in mass fraction, mole fraction or molality</p>
<p style=\"margin-left: 60px;\"><b>refX, refY, refm</b>: Reference composition of species. Needs to be a vector with size of the number of species. The number of species includes the solvent.</p>
<p style=\"margin-left: 60px;\"><br>Additional note: For two phase system there is also the &quot;refLiquidMassFraction&quot; input box. This describes the mass fraction of the mass phase related to the whole mass.</p>
<p><br><b>Tab: Liquid phase</b></p>
<p style=\"margin-left: 30px;\"><b>nL</b>: The number of liquid species including solutes and solvents. Also includes species in their dissociated form. Make sure to include water as solvent.</p>
<p style=\"margin-left: 30px;\"><b>LiquidModel</b>: Choose a liquid model as given in Table 1.</p>
<p style=\"margin-left: 30px;\"><b>substanceNames</b>: A vector containing the names of the substances. Needs to be same size as nL. You can choose any name for each substance. Make sure to include water as solvent (the last element of the vector)</p>
<p style=\"margin-left: 30px;\"><b>data</b>: Data record containing data for calculations like molar Mass, specific gas constant, etc. For further information on data see below.</p>
<p style=\"margin-left: 30px;\"><b>interaction</b>: Binary interaction data record for Gibbs Free Energy Excess models, i.e., Bromley and Pitzer model.</p>
<p><br><b>Tab: Reaction</b></p>
<p style=\"margin-left: 30px;\"><b>nR</b>: Number of gas-liquid and dissociation equilibria.</p>
<p style=\"margin-left: 30px;\"><b>nu</b>: Stoichiometry matrix of gas-liquid and dissociation equilibria. nu has the size of nR x nL.</p>
<h4>Data</h4>
<p>The data for each species can be found in</p>
<p>Gas phase: <a href=\"Media.GasPhase.Common.FluidData\">Media.GasPhase.Common.FluidData</a>, <a href=\"Media.GasPhase.Common.SingleGasesData\">Media.GasPhase.Common.SingleGasesData</a>, <a href=\"Media.GasPhase.Common.MixtureGasesData\">Media.GasPhase.Common.MixtureGasesData</a></p>
<p>Solid phase: <a href=\"Media.SolidPhase.Common.SolidData\">Media.SolidPhase.Common.SolidData</a></p>
<p>Liquid phase: <a href=\"Media.LiquidPhase.Common.SolutesData\">Media.LiquidPhase.Common.SolutesData</a>, <a href=\"Media.LiquidPhase.Common.MixtureSolutesData\">Media.LiquidPhase.Common.MixtureSolutesData</a></p>
<p>There are records for single species. To use data in a mixture of species, a data vector record has to be created as you can see inside the data packages.</p>
</html>"));
end MediumUsage;
