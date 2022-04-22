within ElectrolyteMedia.Media.GasPhase.Common.Functions.IdealGas;
function calc_g_i "Specific Gibbs free energy at temperature T"
  input SI.Temperature T "Temperature";
  output SI.SpecificGibbsFreeEnergy[nGfun] g "Specific Gibbs free energy at temperature T";
protected
  SI.MolarEnergy[nGfun] gm = Molar.calc_g_i(T);
algorithm
  for i in 1:nGfun loop
    g[i] :=gm[i]/datafun[i].MM;
  end for;
//   g :=Molar.calc_g_i(T) ./ datafun.MM;
  annotation (smoothOrder=5);
end calc_g_i;
