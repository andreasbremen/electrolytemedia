within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Bromley.Solvent;
function calc_I_c

  input SI.MassFraction[nLfun] X;

  output Real I_c;

protected
  Real[nLfun-1] mol_i =  calc_mfromX(X);

algorithm
  I_c :=0;
  for i in 1:nLfun-1 loop
    if datafun[i].z > 0 then
      I_c :=I_c + 0.5*mol_i[i]*datafun[i].z^2;
    end if;
  end for;
end calc_I_c;
