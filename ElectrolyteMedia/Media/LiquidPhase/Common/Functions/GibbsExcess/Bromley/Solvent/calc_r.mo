within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Bromley.Solvent;
function calc_r
  "calculates resudue term for mixed electrolyte activity of water"

  input SI.MassFraction[nLfun] X;

  output Real r;
protected
  Real mol_i[nLfun-1] =  calc_mfromX(X);
  Real X_i[nLfun - 1]=calc_X_i(X);
  Real Y_j[nLfun - 1]=calc_Y_j(X);
  Real W_ij[nLfun - 1,nLfun - 1]=calc_W_ij(X);
  Real I=calc_I(X);
  Real I_i[nLfun-1];
algorithm
  for i in 1:nLfun-1 loop
    I_i[i] :=0.5*mol_i[i]*datafun[i].z^2;
  end for;

  r :=0;
  for i in 1:nLfun-1 loop
    for j in 1:nLfun-1 loop
      if datafun[i].z  * datafun[j].z == 0 then
        r :=r;
      else
        r :=r + 0.0156*I*W_ij[i, j]*X_i[i]*Y_j[j]/abs(datafun[i].z)/abs(datafun[j].z);
      end if;
    end for;
    if datafun[i].z == 0 then
      r :=r;
    else
      r :=r - 0.0156*I_i[i]/abs(datafun[i].z);
    end if;
  end for;

end calc_r;
