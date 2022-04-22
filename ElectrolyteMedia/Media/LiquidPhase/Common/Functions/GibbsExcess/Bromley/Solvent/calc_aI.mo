within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Bromley.Solvent;
function calc_aI "calculates 1.5/(z_i*z_j)*I"

  input SI.MassFraction[nLfun] X;

  output Real[nLfun-1,nLfun-1] aI;

protected
  Real I = calc_I(X);

algorithm

  for i in 1:nLfun-1 loop
    for j in 1:nLfun-1 loop
      if datafun[i].z > 0 and datafun[j].z < 0 then
        aI[i,j] := 1.5/abs(datafun[i].z*datafun[j].z)*I;
      end if;
    end for;
  end for;
end calc_aI;
