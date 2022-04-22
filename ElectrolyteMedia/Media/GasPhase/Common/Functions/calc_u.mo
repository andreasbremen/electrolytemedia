within ElectrolyteMedia.Media.GasPhase.Common.Functions;
function calc_u "Calculates specific internal energy of gas phase with corresponding EOS"
  input SI.Temperature T;
  input SI.Density d;
  input SI.MassFraction[:] X;

  output SI.SpecificEnergy u;

protected
  SI.Pressure p;

algorithm
  if GasModelfun == Media.Common.Types.GasModel.Ideal then
    p :=IdealGas.calc_p(T, d, X);
    u :=IdealGas.calc_u(T,p,X);
  elseif GasModelfun == Media.Common.Types.GasModel.PengRobinson then
    u :=PengRobinson.calc_u(T,d,X);
  end if;

end calc_u;
