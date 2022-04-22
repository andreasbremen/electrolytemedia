within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson.Molar;
function calc_b_i
  "Calculates gas minimal volume b_i of each species in Peng Robinson EOS"

  output SI.MolarVolume[nGfun] b_i;
algorithm
   for i in 1:nGfun loop
     b_i[i] :=0.07780*Modelica.Constants.R*datafun[i].T_c/datafun[i].p_c;
   end for;
end calc_b_i;

function calc_B_i
  "Calculates share of gas minimal volume of species in Peng Robinson EOS"

  input SI.MoleFraction[nGfun] y_i;
  output Real[nGfun] B_i(unit="");
protected
  SI.MolarVolume[nGfun] b_i=Molar.calc_b_i();
  SI.MolarVolume b=Molar.calc_b(y_i);
algorithm
  B_i := b_i/(b);
end calc_B_i;
