within ElectrolyteMedia.Media.GasPhase.Common.Functions;
function calc_phi "Calculates fugacity coefficients of gas phase with corresponding EOS"
  input SI.Temperature T;
  input SI.Density d;
  input SI.MassFraction[:] X;

  output Real[nGfun] phi;

algorithm
  if GasModelfun == Media.Common.Types.GasModel.Ideal then
    phi := ones(nGfun);
  elseif GasModelfun == Media.Common.Types.GasModel.PengRobinson then
    phi :=PengRobinson.calc_phi(T,d,X);
  end if;

end calc_phi;
