within ElectrolyteMedia.Media.GasLiquidPhase.Common.MixtureLiquid.Initialization.calc_f;
function GLE_Tp
  "GLE initialization according to Leal et al. (2016): calculation of residuals"
  input Real[nG+nL] x;
  input Real[nX] Xred;
  input SI.Temperature T;
  input SI.Pressure p;
  output Real[nG+nL] f;

protected
  SI.MolarEnergy gi[nG+nL];
  SI.MolarEnergy gim[nG+nL];
  SI.MolarEnergy gr[nR];
  Real[nR] logK;
  Real[nL] gamma_i;
  Real[nL] m;
  Real [nL] a;
  SI.MoleFraction[nG] Yg;
  SI.MoleFraction[nL] Yl;
  SI.MassFraction[nG] Xg;
  SI.MassFraction[nL] Xl;
  Real[nG] phi;
  Real[nG] fug;
  SI.Density dg;
  Real tau = 1e-37;
  SI.MassFraction[nG+nL] z;

  Real[nG+nL] logreacBase;
  Real[nG+nL] x_orig "moles in original order";
algorithm
  assert(sum(x) > 0, "x is zero in GLE_Tp");

  x_orig :=P_to_orig*x;
  z :=tau ./ x_orig;


  Yg :=x_orig[1:nG]/sum(x_orig[1:nG]);
  Yl :=x_orig[1 + nG:nG + nL]/sum(x_orig[1 + nG:nG + nL]);
  Xg :=Functions.calc_X_M(Yg, MMX[1:nG]);
  Xl :=Functions.LiquidFunctions.calc_X(Yl);

  gamma_i := Functions.LiquidFunctions.calc_gamma(T,p,Xl);
  dg :=dg_TpX(T,p,Xg);
  phi :=Functions.GasFunctions.calc_phi(T,dg,Xg);

  m[1:nL-1] :=Functions.LiquidFunctions.calc_mfromY(x_orig[1+ nG:nG + nL]);
  m[nL] :=Xl[nL]/MH2O;

  gi :=Functions.calc_gi(T, p) .* MMX;
  gr :=nu*gi;
  logK :=-gr/Modelica.Constants.R/T;

  a :=gamma_i .* m;
  fug :=Yg .* phi*p/prefg;

  logreacBase :=cat(1,log(fug) - z[1:nG], log(a) - z[nG + 1:nG + nL]);

  //isopotential
  f[1:nR] :=nu*logreacBase - logK;

  //reduced mass balance
  f[nR+1:nF] :=transpose(lambda_id)*x - Xred;

//   //isopotential
//   for r in 1:nR loop
//     f[r] :=(sum({if nu[r,i_] <> 0 then log(x[i_]*phi[i_]*p/prefg)*nu[r,i_] else 0 for i_ in 1:nG})+sum({if nu[r,i_+nG] <> 0 then log(gamma_i[i_]*m[i_])*nu[r, i_+nG] else 0 for i_ in 1:nL}) - log(K[r]));
//   end for;
//
//   //nR+1: closure condition gas phase
//   f[nR+1] :=sum(x[1:nG]) - 1;
//
//   //nR+2: closure condition liquid phase
//   f[nR+2] :=sum(x[1 + nG:nG + nL]) - 1;
//
//   //nR+3: charge balance liquid phase
//   f[nR+3] :=datal[:].z*x[1 + nG:nG + nL - 1];
//
//   //constraints on inert solute molalities and mole fraction of gas phase
//   index :=4;
//    for i in 1:nG+nL loop
//      if sum(abs(nu[:,i])) < Modelica.Constants.eps then
//        if i < nG+1 then
//          f[nR+index] :=1000*(x[i] - inert[i]);
//        else
//          f[nR+index] :=x[i] - x[nG + nL]*MH2O*inert[i];
//        end if;
//        index :=index + 1;
//        if index > nG+nL-nR then
//          break;
//        end if;
//      end if;
//    end for;
end GLE_Tp;
