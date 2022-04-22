within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson.Molar;
function calc_MM "Calculates molar mass of a mixture"

  input SI.MoleFraction y[nGfun];
  output SI.MolarMass MM;

protected
  SI.MolarMass[nGfun] MM_i;

algorithm
  for i in 1:nGfun loop
    MM_i[i] :=datafun[i].MM*y[i];
  end for;
  MM := sum(MM_i);
end calc_MM;
