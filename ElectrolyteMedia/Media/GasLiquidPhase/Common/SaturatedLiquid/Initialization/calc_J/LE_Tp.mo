within ElectrolyteMedia.Media.GasLiquidPhase.Common.SaturatedLiquid.Initialization.calc_J;
function LE_Tp
  "Liquid dissociation equilibrium with constraints on solute molalities"
  input Real[nL] x;
  input SI.Temperature T;
  input SI.Pressure p;
  output Real[nL,nL] J;

protected
  SI.MoleFraction[nL] Y;
  SI.MassFraction[nL] X;
  Real[nL] gamma;
  Real[nL] m;
  SI.MolarMass[nL] MMXl = Functions.LiquidFunctions.calc_MMX();
algorithm
  assert(sum(x)>0,"Input error: x is zero");

  Y :=x/sum(x);

  m[1:nL-1] := Functions.LiquidFunctions.calc_mfromY(Y);
  m[nL] := Y[nL]/MH2O;
  X :=Functions.LiquidFunctions.calc_X(Y);
  gamma := Functions.LiquidFunctions.calc_gamma(T,p,X);

  for i in  1:nL loop
    //isopotential
    for r in 1:nR_L loop
       if i < nL then
         J[r, i] := if abs(nu_L[r,i]) > 0 then nu_L[r,i]/x[i] - nu_L[r,nL]/sum(x) else - nu_L[r,nL]/sum(x);
       else
         J[r,i] :=1/x[i]*(nu_L[r,i]-sum(nu_L[r,1:nL-1])) - nu_L[r,i]/sum(x);
       end if;
    end for;

    //reduced mass balance
      for j in 1:nL - nR_L loop
        J[j + nR_L, i] := lambda_mass_L[i, j]*MMXl[i];
      end for;
  end for;
end LE_Tp;
