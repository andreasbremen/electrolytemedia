within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Bromley.Solvent;
function calc_X_i

  input SI.MassFraction[nLfun] X;

  output Real[nLfun-1] X_i;
protected
  Real[nLfun-1] mol_i =  calc_mfromX(X);
  Real I_i[nLfun-1];
  Real I_c=calc_I_c(X);
algorithm
  for i in 1:nLfun-1 loop
    if datafun[i].z > 0 then
      I_i[i] :=0.5*mol_i[i]*datafun[i].z^2;
    else
      I_i[i] :=0;
    end if;
  end for;

  X_i :=I_i./(I_c+1e-25);

end calc_X_i;
