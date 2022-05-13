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

  diagH[1:ns] :=tau ./ n[1:ns] .^ 2;//z[1:ns]./n[1:ns];//
  diagH[ns+1:ns+nL] :=(1 .- Yl .+ z[ns+1:ns+nL]) ./ n[ns + 1:ns + nL];// + tau ./ n[ns+1:ns+nL] .^ 2;//(1 .- Yl) ./ n[ns + 1:ns + nL] + tau ./ n[ns+1:ns+nL] .^ 2;

  H :=diagonal(diagH);//fill(diagH,ns+nL);//
  H_id :=transpose(P_to_id)*H;

//    for i in 1:ns+nL loop
//      //isopotential
//      for r in 1:nR loop
//        if i < ns+1 then
//           J[r,i] := nu[r,i]*tau/x[i]^2 else 0;//nu[r,i]*z[i]/n[i] else 0;//if abs(nu[r,i]) > 0 then nu[r,i]*z[i]/n[i] else 0;
//         else
//           J[r,i] := if abs(nu[r,i])>0 then nu[r, i]*(1 - Yl[i - ns])/x[i] + nu[r,i]*tau/x[i]^2 else 0;//nu[r, i]*(1 - Yl[i - ns])/n[i] + nu[r,i]*z[i]/n[i] else 0;//if abs(nu[r,i]) > 0 then nu[r, i]*(1 - Yl[i - ns])/n[i] + nu[r,i]*z[i]/n[i] else 0;
//         end if;
//      end for;
//    end for;

  //isopotential
//   for r in 1:nR loop
//     J[r,:] :=nu[r,:].*diagH;
//   end for;
     J[1:nR,:] :=nu*H;//nu_id__ .* H;//nu_id__*H;//
     J[1:nR,:] :=J[1:nR, :]*transpose(P_to_id);//transpose(P_nu_id);

//      J[1:nR,:] :=nu_id_ .* H_id;
//     J[1:nR,:] :=nu_id .* H_id;

  //reduced mass balance
   J[nR+1:nF,:] :=transpose(lambda_id);
end SLE_Tp;
