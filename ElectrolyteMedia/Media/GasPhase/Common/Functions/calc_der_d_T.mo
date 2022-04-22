within ElectrolyteMedia.Media.GasPhase.Common.Functions;
function calc_der_d_T
  "Calculates derivative of density w.r.t. temperature at temperature and pressure with corresponding EOS"
  input SI.Temperature T;
  input SI.Density d;
  input SI.MassFraction[:] X;

  output SI.DerDensityByTemperature ddTp;

protected
  SI.Pressure p;

algorithm
  if GasModelfun == Media.Common.Types.GasModel.Ideal then
    p :=IdealGas.calc_p(T, d, X);
    ddTp :=IdealGas.calc_der_d_i_T(
      T,
      p,
      X);
  elseif GasModelfun == Media.Common.Types.GasModel.PengRobinson then
    ddTp :=PengRobinson.calc_der_d_T(
      T,
      d,
      X);
  end if;

end calc_der_d_T;
