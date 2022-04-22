within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson.Molar;
function calc_a_i
  "Calculates van der Waals attraction of each species in Peng Robinson EOS"
  output Real[nGfun] a_i(unit="N.m4/mol2");
algorithm
  for i in 1:nGfun loop
    a_i[i] :=0.45724*(Modelica.Constants.R*datafun[i].T_c).*(Modelica.Constants.R*datafun[i].T_c)/datafun[i].p_c;
  end for;
end calc_a_i;
