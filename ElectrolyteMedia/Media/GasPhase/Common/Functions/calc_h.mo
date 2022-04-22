within ElectrolyteMedia.Media.GasPhase.Common.Functions;
function calc_h "Calculates specific enthalpy of gas phase with corresponding EOS"
  input SI.Temperature T;
  input SI.Density d;
  input SI.MassFraction[:] X;

  output SI.SpecificEnthalpy h;

protected
  SI.Pressure p;

algorithm
  if GasModelfun == Media.Common.Types.GasModel.Ideal then
    p :=IdealGas.calc_p(T, d, X);
    h :=IdealGas.calc_h(T,p,X);
  elseif GasModelfun == Media.Common.Types.GasModel.PengRobinson then
    h :=PengRobinson.calc_h(T,d,X);
  end if;
annotation(smoothOrder=5);
end calc_h;
