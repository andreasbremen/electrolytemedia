within ElectrolyteMedia.Media.GasPhase.Common.Functions;
function calc_g "Calculates specific gibbs free energy of gas phase with corresponding EOS"
  input SI.Temperature T;
  input SI.Density d;
  input SI.MassFraction[:] X;

  output SI.SpecificEnergy g;

protected
  SI.Pressure p;

algorithm
  if GasModelfun == Media.Common.Types.GasModel.Ideal then
    p :=IdealGas.calc_p(T, d, X);
    g :=IdealGas.calc_g(T,p,X);
  elseif GasModelfun == Media.Common.Types.GasModel.PengRobinson then
    g :=PengRobinson.calc_g(T,d,X);
  end if;

end calc_g;
