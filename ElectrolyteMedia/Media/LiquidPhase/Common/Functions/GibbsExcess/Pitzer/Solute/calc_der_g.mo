within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Solute;
function calc_der_g "Derivative of g function for calculation of der_B_ij"

  input Integer P "1 for alpha1, 2 for alpha2";
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real[nLifun,nLifun] g;
protected
  Real[nLifun,nLifun] x;
  Real[nLifun,nLifun] alpha=if P == 1 then Solute.calc_alpha1() else
      Solute.calc_alpha2();
  Real I = calc_I(X);

algorithm
  if I > 0 then
    if P == 1 then
      for i in 1:nLifun loop
        for j in 1:nLifun loop
          if abs(datafun[i].z) > 0 and abs(datafun[j].z) > 0 then
            x[i,j] :=alpha[i, j]*I^0.5;
            g[i,j] :=-2*(1 - (1 + x[i, j]+x[i,j]^2/2)*exp(-x[i, j]))/x[i, j]^2;
          end if;
        end for;
      end for;
    elseif P == 2 then
      for i in 1:nLifun loop
        for j in 1:nLifun loop
          if abs(datafun[i].z) >= 2 and abs(datafun[j].z) >= 2 then
            x[i,j] :=alpha[i, j]*I^0.5;
            g[i,j] :=-2*(1 - (1 + x[i, j]+x[i,j]^2/2)*exp(-x[i, j]))/x[i, j]^2;
          end if;
        end for;
      end for;
    end if;
  end if;
end calc_der_g;
