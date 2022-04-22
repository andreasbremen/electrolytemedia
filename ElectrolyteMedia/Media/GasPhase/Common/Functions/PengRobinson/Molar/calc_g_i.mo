within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson.Molar;
function calc_g_i "Calculates Gibbs free energy of gas species"

  input SI.Temperature T;
  input SI.Density d;
  input SI.MoleFraction[nGfun] y_i;

  output SI.MolarEnergy g_i[nGfun];
protected
  SI.MolarEnergy[nGfun] g_id=IdealGas.Molar.calc_g_i(T);
  SI.MolarEnergy[nGfun] g_ex= Molar.calc_g_i_ex(T,d,y_i);
algorithm
  g_i :=g_id + g_ex;
end calc_g_i;
