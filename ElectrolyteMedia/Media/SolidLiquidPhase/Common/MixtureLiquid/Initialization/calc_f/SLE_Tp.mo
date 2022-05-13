within ElectrolyteMedia.Media.SolidLiquidPhase.Common.MixtureLiquid.Initialization.calc_f;
function SLE_Tp
  "SLE initialization according to Leal et al. (2016): calculation of residuals"
  input Real[nF] x;
  input Real[nX] Xred;
  input SI.Temperature T;
  input SI.Pressure p;
  output Real[nF] f;

protected
  SI.MolarEnergy gi[ns+nL];
  SI.MolarEnergy gi_id[ns+nL];
  SI.MolarEnergy gr[nR];
  Real[nR] logK;
  Real[nL] gamma_i;
  Real[nL] m;
  Real [nL] a;
  SI.MoleFraction[ns+nL] Y;
  SI.MassFraction[ns+nL] X;
  SI.MoleFraction[nL] Yl;
  SI.MassFraction[nL] Xl;
  SI.MoleFraction[ns] Ys;
  SI.MassFraction[ns] Xs;
  Real tau = 1e-37;
  SI.MassFraction[ns+nL] z;// = {if x[i] > 0 then tau/max(tau,x[i]) else 0 for i in 1:ns+nL};//tau./x;

  Real[ns+nL] logreacBase;
  Real[ns+nL] logreacBase_id;
  Real[ns+nL] x_orig;

  SI.Mass[ns+nL] mass;

algorithm

  //assert(sum(x[1:ns+nL]) > 0, "x is zero in SLE_Tp");

  x_orig :=P_to_orig*x;
  for i in 1:ns+nL loop
    x_orig[i] :=x_orig[i];//max(tau, x_orig[i]);//x_orig without pivoting, max formulation with pivoting
  end for;
//     z :={if x_orig[i] > 0 then tau/max(1e-20*tau, x_orig[i]) else 0 for i in 1:ns + nL};
       z :=tau ./ x_orig;

  Y :=x_orig[1:ns+nL]/sum(x_orig[1:ns+nL]);
  X :=Functions.calc_Xfull(Y);

  Ys :=Y[1:ns]/sum(Y[1:ns]);
  Yl :=Y[1 + ns:ns + nL]/sum(Y[1 + ns:ns + nL]);

  Xs :=Functions.calc_X_M(Ys, MMX[1:ns]);
  Xl :=Functions.calc_X_M(Yl, MMX[1 + ns:ns + nL]);

  m[1:nL-1] := Functions.LiquidFunctions.calc_mfromY(Yl);
  m[nL] := Yl[nL]/MH2O;

  gamma_i := Functions.LiquidFunctions.calc_gamma(T,p,Xl);
  gi := Functions.calc_gi(T, p) .* MMX;
  gi_id :=P_to_id*gi;

//    gr :=nu_id*gi_id;
   gr :=nu*gi;
  logK :=-gr/Modelica.Constants.R/T;

  a :=gamma_i .* m;
  logreacBase :=cat(1,-z[1:ns],log(a)-z[ns+1:ns+nL]);
   logreacBase_id :=P_to_id*logreacBase;

//   x_id :=P_to_id*x;

  //isopotential
   f[1:nR] :=nu*logreacBase - logK;
//   f[1:nR] :=nu_id*logreacBase_id - logK;

  //reduced mass balance
  // f[nR+1:nF] :=transpose(lambda_mass)*mass - Xred;
//   f[nR+1:nF] :=transpose(lambda_nu)*x_id - Xred;//transpose(lambda_mass)*mass - Xred;
  f[nR+1:nF] :=transpose(lambda_id)*x - Xred;//transpose(lambda_mass)*mass - Xred;
end SLE_Tp;
