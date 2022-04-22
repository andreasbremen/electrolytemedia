within ElectrolyteMedia.Media.GasLiquidPhase.Common.SaturatedLiquid;
function X_vfracTXred
  "Liquid dissociation equlibrium at bubble point temperature T, with vapor fraction vfrac and reduced mass fractions (sum(Xred)=1) as input. Output is quasi-homogeneous mass fraction of gas and liquid species."
   input Real vfrac;
   input SI.Temperature T;
   input MassFraction[nX] Xred;
   input MassFraction[nF] Xinit = zeros(nF);
   output Real[nG+nL] X;

protected
  Real[nX] Xredl;
  Real eps = 1e-6;
  Boolean solutionfound = false;
  SI.Pressure p;
  SI.MassFraction[nG] Xg;
  SI.MassFraction[nL] Xl1;
  SI.MassFraction[nL] Xl2;

  parameter Real init = 1e-12;
  parameter Integer[nR] firstnonzero = {Modelica.Math.BooleanVectors.firstTrueIndex({nu[r,i] <> 0 for i in 1:nF}) for r in 1:nR};
  Real[nR,nF] B;
  Real[nR] b;
  SI.MassFraction[nF] Xfull_in;
algorithm

  // find initial guess Xfull_in with positive elements only
  for r in 1:nR loop
    B[r,firstnonzero[r]] :=1;
    b[r] :=0;
  end for;
  Xfull_in :=Modelica.Math.Matrices.equalityLeastSquares(transpose(lambda_mass), Xred,B,b);
  Xfull_in :=Xfull_in/sum(Xfull_in);

  // set vapor fraction
  Xfull_in[nG+nL] :=Xfull_in[nG + nL] + Xfull_in[1] - vfrac;
  Xfull_in[1] :=vfrac;

  Xg :=fill(1/nG, nG);

  Xredl :=transpose(lambda_mass_L)*Xfull_in[nG+1:nF]/sum(transpose(lambda_mass_L)*Xfull_in[nG+1:nF]);

  // initial guess for Newton solver
  if sum(Xinit)> 0 then
    Xl1 :=Xinit[nG+1:nF]/sum(Xinit[nG+1:nF]);
  else
    Xl1 :=Xfull_in[nG+1:nF]/sum(Xfull_in[nG+1:nF]);
    for i in 1:nL loop
      Xl1[i] :=max(init, Xfull_in[nG+i]);
    end for;
  end if;

//   Xl1 :=Xfull_in[nG+1:nF]/sum(Xfull_in[nG+1:nF]);
//   Xredl :=transpose(lambda_mass_L)*Xl1/sum(transpose(lambda_mass_L)*Xl1);
//   for i in 1:nL loop
//     if Xl1[i] < init then
//       Xl1[i] :=init;
//     end if;
//   end for;

  X :=cat(1,vfrac*Xg,(1 - vfrac)*Xl1);

  while not solutionfound loop
    p :=p_TX(T, X);
    Xl2 :=X_pTXred(p, T, Xredl, Xl1);
    if Modelica.Math.Vectors.norm((Xl1-Xl2)./Xl1,Modelica.Constants.inf)< eps then
      solutionfound :=true;
    end if;
    Xl1 :=Xl2;
    X :=cat(1,vfrac*Xg,(1 - vfrac)*Xl1);
  end while;

end X_vfracTXred;
