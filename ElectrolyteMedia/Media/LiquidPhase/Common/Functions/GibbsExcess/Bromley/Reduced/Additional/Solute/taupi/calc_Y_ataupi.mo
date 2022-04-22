within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Bromley.Reduced.Additional.Solute.taupi;
function calc_Y_ataupi "Helper function to calculate Gibbs derivative"

  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nLfun] X;

  output Real[nLfun-1] Y_ataupi;
protected
  Real[nLfun - 1,nLfun - 1] loggamma_ik_0taupi=calc_loggamma_ik_0taupi(T,p,X);
  Real I=calc_I(X);
  Real[nLfun-1] mol_i =  calc_mfromX(X);
algorithm
  Y_ataupi :=zeros(nLfun-1);
  for i in 1:nLfun-1 loop
    if datafun[i].z < 0 then
      for k in 1:nLfun-1 loop
        if datafun[k].z > 0 then
          Y_ataupi[i] :=Y_ataupi[i] + (abs(datafun[i].z) + abs(datafun[k].z))^2/4*mol_i[k]/I * loggamma_ik_0taupi[i,k];
        end if;
      end for;
    end if;
  end for;
end calc_Y_ataupi;
