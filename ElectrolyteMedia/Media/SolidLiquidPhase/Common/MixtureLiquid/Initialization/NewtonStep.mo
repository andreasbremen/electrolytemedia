within ElectrolyteMedia.Media.SolidLiquidPhase.Common.MixtureLiquid.Initialization;
function NewtonStep
  "Calculates Newton step with step size such that x[t,i] > 0"
  input Real[nNS] x;
  input Real[nNS] f;
  input Real[nNS,nNS] J;
  input Integer nNS;
  output Real[nNS] x_;

protected
  Real beta=1;
  Real[nNS] alpha;
  Real delta = 0.9999;
   Real[nNS,nNS] invJ= Modelica.Math.Matrices.inv(J);
  Real[nNS] delta_x =  invJ*f;//Modelica.Math.Matrices.solve(J,f);//
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

//    while min(x_temp) < 0 loop
//      x_temp :=x - 1/beta*invJ*f;
//      beta :=beta *10;
//    end while;
//    assert(min(x_temp)>0,"Newton step too small");
//    x_ :=x_temp;

end NewtonStep;
