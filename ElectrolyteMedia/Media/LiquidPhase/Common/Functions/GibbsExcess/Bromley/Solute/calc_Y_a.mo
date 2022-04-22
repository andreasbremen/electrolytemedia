within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Bromley.Solute;
function calc_Y_a "Term for calculation of loggamma"

  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nLfun] X;

  output Real[nLfun-1] Y_a;
protected
  Real[nLfun - 1,nLfun - 1] loggamma_ik_0=calc_loggamma_ik_0( T,p,X);
  Real I=calc_I(X);
  Real[nLfun-1] mol_i = calc_mfromX(X);
algorithm
  Y_a :=zeros(nLfun-1);
  if I > 0 then
    for i in 1:nLfun-1 loop
      if datafun[i].z < 0 then
        for k in 1:nLfun-1 loop
          if datafun[k].z > 0 then
            Y_a[i] :=Y_a[i] + (abs(datafun[i].z) + abs(datafun[k].z))^2/4*mol_i[k]/I * loggamma_ik_0[i,k];
          end if;
        end for;
      end if;
    end for;
  end if;
end calc_Y_a;
