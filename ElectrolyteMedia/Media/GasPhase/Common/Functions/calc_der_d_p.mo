within ElectrolyteMedia.Media.GasPhase.Common.Functions;
function calc_der_d_p
  "Calculates derivative of density w.r.t. pressure at temperature and pressure with corresponding EOS"
  input SI.Temperature T;
  input SI.Density d;
  input SI.MassFraction[nGfun] X;

  output SI.DerDensityByPressure ddpT;

protected
  SI.Pressure p;

algorithm
  if GasModelfun == Media.Common.Types.GasModel.Ideal then
    p :=IdealGas.calc_p(T, d, X);
    ddpT :=IdealGas.calc_der_d_i_p(
      T,
      p,
      X);
  elseif GasModelfun == Media.Common.Types.GasModel.PengRobinson then
    ddpT :=PengRobinson.calc_der_d_p(
      T,
      d,
      X);
  end if;

end calc_der_d_p;
