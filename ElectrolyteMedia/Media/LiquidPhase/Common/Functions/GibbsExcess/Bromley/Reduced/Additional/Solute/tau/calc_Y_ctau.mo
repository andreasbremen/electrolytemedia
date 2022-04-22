within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Bromley.Reduced.Additional.Solute.tau;
function calc_Y_ctau "Helper function to calculate Gibbs derivative"

  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nLfun] X;

  output Real[nLfun-1] Y_ctau;
protected
  Real[nLfun - 1,nLfun - 1] loggamma_ik_0tau=calc_loggamma_ik_0tau(T,p,X);
  Real I=calc_I(X);
  Real[nLfun-1] mol_i = calc_mfromX(X);
algorithm
  Y_ctau :=zeros(nLfun-1);
  for i in 1:nLfun-1 loop
    if datafun[i].z > 0 then
      for k in 1:nLfun-1 loop
        if datafun[k].z < 0 then
          Y_ctau[i] :=Y_ctau[i] + (abs(datafun[i].z) + abs(datafun[k].z))^2/4*mol_i[k]/I * loggamma_ik_0tau[i,k];
        end if;
      end for;
    end if;
  end for;
end calc_Y_ctau;
