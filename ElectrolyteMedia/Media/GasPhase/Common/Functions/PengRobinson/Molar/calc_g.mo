within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson.Molar;
function calc_g
  "Calculates molar Gibbs free energy of gas phase with Peng Robinson"
  input SI.Temperature T;
  input SI.Density d;
  input SI.MoleFraction[nGfun] y_i;

  output SI.MolarEnergy g;
protected
  SI.MolarEnergy[nGfun] g_i=Molar.calc_g_i(T,d,y_i);
algorithm
  g := g_i*y_i;
end calc_g;
