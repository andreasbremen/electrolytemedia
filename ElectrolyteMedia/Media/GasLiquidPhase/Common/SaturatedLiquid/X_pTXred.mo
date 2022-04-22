within ElectrolyteMedia.Media.GasLiquidPhase.Common.SaturatedLiquid;
function X_pTXred
  "Liquid dissociation equilibrium calculated via Newton method with homotopy continuation method"
  input SI.Pressure p;
  input SI.Temperature T;
  input MassFraction[nX_L] Xredl;
  input MassFraction[nL] Xinit = zeros(nL);
  output Real[nL] Xfull;

protected
  parameter Real init = 1e-10;
  parameter Real eps = 1e-5;
  parameter Integer[nR_L] firstnonzero = {Modelica.Math.BooleanVectors.firstTrueIndex({nu_L[r,i] <> 0 for i in 1:nL}) for r in 1:nR_L};
  parameter Integer[nR_L] secondnonzero = {Modelica.Math.BooleanVectors.firstTrueIndex({if i>firstnonzero[r] then nu_L[r,i] <> 0 else false for i in 1:nL}) for r in 1:nR_L};

  Real[nL,nL] J;
  Real[nL] f;
  Real[nL] x1;
  Real[nL] x2;

  Real[nR_L,nL] B;
  Real[nR_L] b;

  SI.MassFraction[nL] Xfull_in;
  SI.Mass[nL] mfull_in;
  SI.MoleFraction[nL] Yfull_in;
  SI.AmountOfSubstance[nL] nfull_in;

  SI.MoleFraction[nL] Yfull;

  SI.MoleFraction[nL] Yinit;
  SI.AmountOfSubstance[nL] ninit;
  SI.Mass[nL] minit;
  Real scale;

  Boolean solutionfound = false;

  SI.MolarMass[nL] MMXl = Functions.LiquidFunctions.calc_MMX();

algorithm
  // find initial guess x1 for Newton solver
   for r in 1:nR_L loop
     B[r,firstnonzero[r]] :=1;
     b[r] :=0;
   end for;
   Xfull_in :=Media.Common.equalityLeastSquares(transpose(lambda_mass_L),Xredl,B,b);
   if sum(Xinit)> 0 then
   // ninit from initial guess provided as input (Xinit)
     scale :=sum(transpose(lambda_mass_L)*Xinit);
     minit :=Xinit*scale;
     Yinit :=Functions.LiquidFunctions.calc_Y(Xinit);
     ninit :=minit./MMXl;
     for i in 1:nL loop
       x1[i] :=max(init, ninit[i]);
     end for;
   elseif sum(Xfull_in) > 0 then
   // nfull_in with positive elements only to fulfill transpose(lambda_mass)*mfull_in = Xred;
     mfull_in :=Xfull_in;
     nfull_in :=mfull_in ./ MMXl;
     for i in 1:nL loop
       x1[i] :=max(init, nfull_in[i]);
     end for;
   else
   //else generic initial guess
     x1 :=fill(1/nL, nL);
   end if;

  //Newton algorithm
  solutionfound :=false;
  while not solutionfound loop
    f :=Initialization.calc_f.LE_Tp(x1,Xredl,T,p);
    J :=Initialization.calc_J.LE_Tp(x1,T,p);
    x2 := Initialization.NewtonStep(x1,f,J,nL);
    if Modelica.Math.Vectors.norm(f,Modelica.Constants.inf) < eps then
      solutionfound :=true;
    end if;
    x1 :=x2;
  end while;

  Yfull :=x2/sum(x2);
  Xfull :=Functions.LiquidFunctions.calc_X(Yfull);

//   // find initial guess nfull_in with positive elements only and 1kg in total
//   for r in 1:nR_L loop
//     B[r,firstnonzero[r]] :=1;
//     b[r] :=0;
//   end for;
//   Xfull_in :=Modelica.Math.Matrices.equalityLeastSquares(transpose(lambda_mass_L), Xred,B,b);
//   Xfull_in :=Xfull_in/sum(Xfull_in);
//   Yfull_in :=Functions.LiquidFunctions.calc_Y(Xfull_in);
//   nfull_in :=Yfull_in/MH2O;
//
//   nred :=transpose(lambda_L)*nfull_in;
//
//   // initial guess for Newton solver
//   if sum(Xinit)> 0 then
//     Yinit :=Functions.LiquidFunctions.calc_Y(Xinit);
//     ninit :=Yinit/MH2O;
//     x1 :=ninit;
//   else
//     for i in 1:nL loop
//       x1[i] :=max(init, nfull_in[i]);
//     end for;
//   end if;
//
//   // Newton algorithm
//   solutionfound :=false;
//   while not solutionfound loop
//     f :=Initialization.calc_f.LE_Tp(x1,nred,T,p);
//     J :=Initialization.calc_J.LE_Tp(x1,T,p);
//     x2 := Initialization.NewtonStep(x1,f,Modelica.Math.Matrices.inv(J),nL);
//     if Modelica.Math.Vectors.norm(f,Modelica.Constants.inf) < eps then
//       solutionfound :=true;
//     end if;
//     x1 :=x2;
//   end while;
//
//   Yfull :=x2/sum(x2);
//   Xfull :=Functions.LiquidFunctions.calc_X(Yfull);
end X_pTXred;
