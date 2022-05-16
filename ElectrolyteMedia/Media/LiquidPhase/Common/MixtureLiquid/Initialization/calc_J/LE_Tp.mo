within ElectrolyteMedia.Media.LiquidPhase.Common.MixtureLiquid.Initialization.calc_J;
function LE_Tp
  "Liquid dissociation equilibrium with constraints on solute molalities"
  input Real[nL] x;
  input SI.Temperature T;
  input SI.Pressure p;
  output Real[nL,nL] J;

protected
  SI.AmountOfSubstance[nF] n;
  SI.MoleFraction[nL] Y;
  SI.MassFraction[nL] X;
  Real[nL] gamma_i;
  Real[nL] m;
  Real tau = 1e-37;
  SI.MassFraction[nF] z;
  Real[nF,nF] H;
  Real[nF,nF] H_id;
  Real[nF] diagH;
  Real[nF] x_orig;
algorithm

  x_orig :=P_to_orig*x;

  for i in 1:nF loop
    n[i] :=x_orig[i];//max(tau, x_orig[i]);//x_orig without pivoting, max formulation with pivoting
  end for;

  Y :=n[1:nF]/sum(n[1:nF]);
  X :=Functions.calc_X(Y);

  z :=tau./n;

  diagH :=(1 .- Y .+ z) ./ n;

  H :=diagonal(diagH);
  H_id :=transpose(P_to_id)*H;

  //isopotential
  J[1:nR,:] :=nu*H;
  J[1:nR,:] :=J[1:nR, :]*transpose(P_to_id);

  //reduced mass balance
  J[nR+1:nF,:] :=transpose(lambda_id);
end LE_Tp;
