within ElectrolyteMedia.Media.SolidPhase.Common.Functions.Molar;
function calc_xi "Calculates Einstein function"

  input SI.Temperature T;

  output Real xi[nSfun];
protected
  Real u[nSfun]=Molar.calc_theta_T(T);
algorithm
  xi := u .^ 2. * exp(u) ./ (exp(u) - ones(nSfun)) .^ 2;
  annotation(smoothOrder=20);
end calc_xi;
