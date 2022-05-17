within ElectrolyteMedia.Media.GasLiquidPhase.Common.MixtureLiquid.Initialization.calc_J;
function GLE_Tp
  "SLE initialization according to Leal et al. (2016): calculation of Jacobian"
  input Real[nF] x;
  input SI.Temperature T;
  input SI.Pressure p;
  output Real[nF,nF] J;

protected
  SI.AmountOfSubstance[nG+nL] n;
  SI.MoleFraction[nG+nL] Y;
  SI.MassFraction[nG+nL] X;
  SI.MoleFraction[nL] Yl;
  SI.MassFraction[nL] Xl;
  SI.MoleFraction[nG] Yg;
  SI.MassFraction[nG] Xg;
  Real tau = 1e-37;
  SI.MassFraction[nG+nL] z;
  Real[nG+nL,nG+nL] H;
  Real[nG+nL,nG+nL] H_id;
  Real[nG+nL] diagH;
  Real[nG+nL] x_orig;
algorithm

  x_orig :=P_to_orig*x;

  for i in 1:nG+nL loop
    n[i] :=x_orig[i];//max(tau, x_orig[i]);//x_orig without pivoting, max formulation with pivoting
  end for;

  Y :=n[1:nG+nL]/sum(n[1:nG+nL]);
  X :=Functions.calc_Xfull(Y);

  Yg :=Y[1:nG]/sum(Y[1:nG]);
  Yl :=Y[1 + nG:nG + nL]/sum(Y[1 + nG:nG + nL]);

  z :=tau./n;

  diagH[1:nG] :=(1 .- Yg .+ z[1:nG]) ./ n[1:nG];
  diagH[nG+1:nG+nL] :=(1 .- Yl .+ z[nG+1:nG+nL]) ./ n[nG + 1:nG + nL];

  H :=diagonal(diagH);
  H_id :=transpose(P_to_id)*H;

  //isopotential
  J[1:nR,:] :=nu*H;
  J[1:nR,:] :=J[1:nR, :]*transpose(P_to_id);

  //reduced mass balance
  J[nR+1:nF,:] :=transpose(lambda_id);
end GLE_Tp;
