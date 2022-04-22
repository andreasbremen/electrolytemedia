within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Bromley.Solute;
function calc_X_c "Term for calculation of loggamma"

  input SI.MassFraction[nLfun] X;

  output Real[nLfun-1] X_c;

protected
  Real I=calc_I(X);
  Real[nLfun-1] mol_i = calc_mfromX(X);

algorithm
  X_c :=zeros(nLfun-1);
  if I > 0 then
    for i in 1:nLfun-1 loop
      if datafun[i].z > 0 then
        for k in 1:nLfun-1 loop
          if datafun[k].z < 0 then
            X_c[i] :=X_c[i] + abs(datafun[i].z*datafun[k].z)*(abs(datafun[i].z) + abs(datafun[k].z))^2/4*mol_i[k]/I;
          end if;
        end for;
      end if;
    end for;
  end if;
end calc_X_c;
