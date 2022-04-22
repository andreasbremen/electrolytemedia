within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Bromley.Solute;
function calc_X_a "Term for calculation of loggamma"

  input SI.MassFraction[nLfun] X;

  output Real[nLfun-1] X_a;

protected
  Real I=calc_I(X);
  Real[nLfun-1] mol_i = calc_mfromX(X);
algorithm
  X_a :=zeros(nLfun-1);
  if I > 0 then
    for j in 1:nLfun-1 loop
      if datafun[j].z < 0 then
        for k in 1:nLfun-1 loop
          if datafun[k].z > 0 then
            X_a[j] :=X_a[j] + abs(datafun[j].z*datafun[k].z)*(abs(datafun[j].z) + abs(datafun[k].z))^2/4*mol_i[k]/I;
          end if;
        end for;
      end if;
    end for;
  end if;
end calc_X_a;
