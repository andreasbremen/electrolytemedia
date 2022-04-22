within ElectrolyteMedia.Media.GasLiquidPhase.Common.MixtureLiquid.Initialization;
function NewtonStep
  "Calculates Newton step with step size such that x_ > 0"
  input Real[nNS] x;
  input Real[nNS] f;
  input Real[nNS,nNS] J;
  input Integer nNS;
  output Real[nNS] x_;
protected
  Real[nNS] x_temp = fill(-1,nNS);
  Real beta=1;
  Real[nNS] alpha;
  Real delta = 0.9999;
  Real[nNS,nNS] invJ= Modelica.Math.Matrices.inv(J);
  Real[nNS] delta_x = invJ*f;
algorithm
  for i in 1:nNS loop
    if x[i] - delta_x[i] < 0 then
      alpha[i] :=delta*x[i]/delta_x[i];
    else
      alpha[i] :=1;
    end if;
  end for;
  beta :=min(alpha);
  x_ :=x - beta*delta_x;
end NewtonStep;
