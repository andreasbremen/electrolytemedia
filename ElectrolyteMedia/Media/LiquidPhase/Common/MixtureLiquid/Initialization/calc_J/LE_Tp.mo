within ElectrolyteMedia.Media.LiquidPhase.Common.MixtureLiquid.Initialization.calc_J;
function LE_Tp
  "Liquid dissociation equilibrium with constraints on solute molalities"
  input Real[nL] x;
  input SI.Temperature T;
  input SI.Pressure p;
  output Real[nL,nL] J;

protected
  SI.MoleFraction[nL] Y;
  SI.MassFraction[nL] X;
  Real[nL] gamma_i;
  Real[nL] m;
algorithm
  assert(sum(x)>0,"Input error: x is zero");
  Y :=x/sum(x);

  m[1:nL-1] := Functions.calc_mfromY(Y);
  m[nL] := Y[nL]/MH2O;

  X :=Functions.calc_X(Y);

  gamma_i := Functions.calc_gamma(T,p,X);

  for i in  1:nL loop
    //isopotential
    for r in 1:nR loop
      if i < nL then
        J[r, i] := if abs(nu[r,i]) > 0 then nu[r,i]/x[i] - nu[r,nL]/sum(x) else - nu[r,nL]/sum(x);
      else
        J[r,i] :=1/x[i]*(nu[r,i]-sum(nu[r,1:nL-1])) - nu[r,i]/sum(x);
      end if;
    end for;

    //reduced mass balance
    for j in 1:nL - nR loop
      J[j + nR, i] := lambda_mass[i, j]*MMX[i];
    end for;
  end for;
end LE_Tp;
