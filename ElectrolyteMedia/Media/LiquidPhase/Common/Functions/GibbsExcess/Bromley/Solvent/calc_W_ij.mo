within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Bromley.Solvent;
function calc_W_ij

  input SI.MassFraction[nLfun] X;

  output Real[nLfun-1,nLfun-1] W_ij;
protected
  Real X_i[nLfun - 1]=calc_X_i(X);
  Real Y_j[nLfun - 1]=calc_Y_j(X);
  Real I_c=calc_I_c(X);
  Real I_a=calc_I_a(X);
  Real I=calc_I(X);
algorithm
  for i in 1:nLfun-1 loop
    for j in 1:nLfun-1 loop
      if datafun[i].z > 0 and datafun[j].z < 0 then
        W_ij[i,j] :=X_i[i]*Y_j[j]*(abs(datafun[i].z) + abs(datafun[j].z))^2/abs(datafun[i].z*datafun[j].z)
          *I_c*I_a/(I^2+1e-25);
      else
        W_ij[i,j] :=0;
      end if;
    end for;
  end for;

end calc_W_ij;
