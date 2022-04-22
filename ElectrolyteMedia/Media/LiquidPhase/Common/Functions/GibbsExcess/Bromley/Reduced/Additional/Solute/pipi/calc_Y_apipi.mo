within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Bromley.Reduced.Additional.Solute.pipi;
function calc_Y_apipi "Helper function to calculate Gibbs derivative"

  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nLfun] X;

  output Real[nLfun-1] Y_apipi;
protected
  Real[nLfun - 1,nLfun - 1] loggamma_ik_0pipi=calc_loggamma_ik_0pipi(T,p,X);
  Real I=calc_I(X);
  Real[nLfun-1] mol_i =  calc_mfromX(X);
algorithm
  Y_apipi :=zeros(nLfun-1);
  for i in 1:nLfun-1 loop
    if datafun[i].z < 0 then
      for k in 1:nLfun-1 loop
        if datafun[k].z > 0 then
          Y_apipi[i] :=Y_apipi[i] + (abs(datafun[i].z) + abs(datafun[k].z))^2/4*mol_i[k]/I * loggamma_ik_0pipi[i,k];
        end if;
      end for;
    end if;
  end for;
end calc_Y_apipi;
