within ElectrolyteMedia.Media.SolidLiquidPhase.Common.MixtureLiquid.Initialization.calc_J;
function SLE_Tp
  "SLE initialization according to Leal et al. (2016): calculation of Jacobian"
  input Real[nF] x;
  input SI.Temperature T;
  input SI.Pressure p;
  output Real[nF,nF] J;

protected
  SI.AmountOfSubstance[ns+nL] n;
  SI.MoleFraction[ns+nL] Y;
  SI.MassFraction[ns+nL] X;
  SI.MoleFraction[nL] Yl;
  SI.MassFraction[nL] Xl;
  SI.MoleFraction[ns] Ys;
  SI.MassFraction[ns] Xs;
  Real tau = 1e-37;
  SI.MassFraction[ns+nL] z;
  Real[ns+nL,ns+nL] H;
  Real[ns+nL,ns+nL] H_id;
  Real[ns+nL] diagH;
  Real[ns+nL] x_orig;
algorithm

  x_orig :=P_to_orig*x;

  for i in 1:ns+nL loop
    n[i] :=x_orig[i];//max(tau, x_orig[i]);//x_orig without pivoting, max formulation with pivoting
  end for;

  Y :=n[1:ns+nL]/sum(n[1:ns+nL]);
  X :=Functions.calc_Xfull(Y);

  Ys :=Y[1:ns]/sum(Y[1:ns]);
  Yl :=Y[1 + ns:ns + nL]/sum(Y[1 + ns:ns + nL]);

  z :=tau./n;

  diagH[1:ns] :=tau ./ n[1:ns] .^ 2;
  diagH[ns+1:ns+nL] :=(1 .- Yl .+ z[ns+1:ns+nL]) ./ n[ns + 1:ns + nL];

  H :=diagonal(diagH);
  H_id :=transpose(P_to_id)*H;

  //isopotential
  J[1:nR,:] :=nu*H;
  J[1:nR,:] :=J[1:nR, :]*transpose(P_to_id);

  //reduced mass balance
  J[nR+1:nF,:] :=transpose(lambda_id);
end SLE_Tp;
