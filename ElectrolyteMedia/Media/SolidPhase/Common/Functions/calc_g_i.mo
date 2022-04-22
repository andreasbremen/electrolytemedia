within ElectrolyteMedia.Media.SolidPhase.Common.Functions;
function calc_g_i "Calculates Gibbs free energy of mineral species"

  input SI.Temperature T;
  input SI.Pressure p;

  output SI.SpecificEnergy g[nSfun];

protected
  SI.MolarEnergy[nSfun] g_m=Molar.calc_g_i(T, p);

algorithm

  for i in 1:nSfun loop
    g[i] :=g_m[i]/datafun[i].MM;
  end for;

  annotation(smoothOrder=5);
end calc_g_i;
