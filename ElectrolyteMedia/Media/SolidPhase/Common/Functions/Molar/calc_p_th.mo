within ElectrolyteMedia.Media.SolidPhase.Common.Functions.Molar;
function calc_p_th
  "Calculates p_th for mineral EOS from Holland & Powell 2011"

  input SI.Temperature T;

  output SI.Pressure p_th[nSfun];
protected
  parameter Real xi_0[nSfun]=Molar.calc_xi(T_0);
  Real u[nSfun]=Molar.calc_theta_T( T);
  parameter Real u_0[nSfun]=Molar.calc_theta_T(T_0);
  parameter Real theta_0[nSfun]=Molar.calc_theta_0();
algorithm
  p_th := datafun[:].alpha_0 .* datafun[:].k_0 .* theta_0 ./ xi_0 .* (ones(nSfun) ./ (exp(u) - ones(nSfun)) - ones(nSfun) ./ (exp(u_0) - ones(nSfun)));
  annotation(smoothOrder=5);
end calc_p_th;
