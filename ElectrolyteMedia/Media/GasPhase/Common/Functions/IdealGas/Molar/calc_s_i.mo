within ElectrolyteMedia.Media.GasPhase.Common.Functions.IdealGas.Molar;
function calc_s_i
  "Calculates molar enthalpy of ideal gas phase at reference pressure"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.MoleFraction[nGfun] y;
  output Modelica.SIunits.MolarEntropy s;
  //Modelica.SIunits.MoleFraction[nGfun] y = cat(1,y_i,{1-sum(y_i)});
protected
  Modelica.SIunits.MolarEntropy s_i[nGfun]=calc_s0_i(T);
algorithm

  s :=s_i*y- y .* Modelica.Constants.R*log(y);
  annotation (Inline=false,smoothOrder=2);
end calc_s_i;
