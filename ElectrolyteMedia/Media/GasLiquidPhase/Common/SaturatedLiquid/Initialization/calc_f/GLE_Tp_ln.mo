within ElectrolyteMedia.Media.GasLiquidPhase.Common.SaturatedLiquid.Initialization.calc_f;
function GLE_Tp_ln
  "Gas-liquid and dissociation equilibrium with constraints on solute molalities and moles of gas phase"
  input Real[nG+nL] x "mole fractions of gas and liquid phase for each phase separately";
  input Real[nG+nL] inert "mole fraction of inert gas phase species and molalities of inert liquid phase species";
  input SI.Temperature T;
//   input Real vfrac;
  output Real[nG+nL] f;

protected
  SI.MoleFraction[nG] Yg = x[1:nG]/sum(x[1:nG]);
  SI.MoleFraction[nL] Yl = x[1+nG:nG+nL]/sum(x[1+nG:nG+nL]);
  SI.MassFraction[nG] Xg = Functions.calc_X_M(Yg, MMX[1:nG]);
  SI.MassFraction[nL] Xl = Functions.LiquidFunctions.calc_X(Yl);
  SI.Pressure p=p_TX(T, cat(
      1,
      Xg,
      Xl));
  SI.SpecificEnergy gi[nG+nL];
  SI.MolarEnergy gim[nG+nL];
  SI.MolarEnergy gr[nR];
  Real[nR] K;
  Real[nL] gamma_i;
  Real[nL] m;
  Real[nG] phi;
  Integer index = 1;
  SI.Density dg;
algorithm
  assert(sum(x) > 0, "x is zero in GLE_Tp");
  gi :=Functions.calc_gi(T, p);
  for i in 1:nG+nL loop
    gim[i] :=gi[i]*MMX[i];
  end for;

  gamma_i := Functions.LiquidFunctions.calc_gamma(T,p,Xl);                                                  // cat(1,ones(nL-1),{MH2O});//cat(1,ones(nL-1),{MH2O});//ForsteriteCarbonationModel.Functions.Aqueous.calc_gamma(T,rho,x[1:nG+nL - 1],params.datal_aqs);
  dg :=dg_TpX(T,p,Xg);
  phi :=Functions.GasFunctions.calc_phi(T,dg,Xg);//fill(1, nG);//

  m[1:nL-1] :=Functions.LiquidFunctions.calc_mfromY(x[1+ nG:nG + nL]);
  m[nL] :=x[nG+nL]/MH2O;

  for j in 1:nR loop
    gr[j] :=sum({nu[j, i]*gim[i] for i in 1:nG+nL});
    K[j] :=exp(-gr[j]/Modelica.Constants.R/T);
  end for;

  //isopotential
  for r in 1:nR loop
    f[r] :=(sum({if nu[r,i_] <> 0 then log(x[i_]*phi[i_]*p/prefg)*nu[r,i_] else 0 for i_ in 1:nG})+sum({if nu[r,i_+nG] <> 0 then log(gamma_i[i_]*m[i_])*nu[r, i_+nG] else 0 for i_ in 1:nL}) - log(K[r]));
  end for;

  //nR+1: closure condition gas phase
//   f[nR+1] :=sum(x[1:nG]) - 1;

  //nR+2: closure condition liquid phase
  f[nR+1] :=sum(x[1 + nG:nG + nL]) - 1;

  //nR+3: charge balance liquid phase
  f[nR+2] :=datal[:].z*x[1 + nG:nG + nL - 1];

  //constraints on inert solute molalities and mole fraction of gas phase
  index :=3;
   for i in 1:nG+nL loop
     if sum(abs(nu[:,i])) < Modelica.Constants.eps then
       if i < nG+1 then
         f[nR+index] :=1000*(x[i] - inert[i]);
       else
         f[nR+index] :=x[i] - x[nG + nL]*MH2O*inert[i];
       end if;
       index :=index + 1;
       if index > nG+nL-nR then
         break;
       end if;
     end if;
   end for;
end GLE_Tp_ln;
