within ElectrolyteMedia.Media.LiquidPhase.Common.MixtureLiquid;
function X_pTXred
  "Liquid dissociation equilibrium calculated via Newton method with homotopy continuation method"
  input SI.Pressure p;
  input SI.Temperature T;
  input MassFraction[nX] Xred;
  input MassFraction[nL] Xinit = zeros(nL);
  output Real[nL] Xfull;

protected
  parameter Real init = 1e-10;
  parameter Real eps = 1e-6;
  parameter Integer[nR] firstnonzero = {Modelica.Math.BooleanVectors.firstTrueIndex({nu[r,i] <> 0 for i in 1:nL}) for r in 1:nR};
  parameter Integer[nR] secondnonzero = {Modelica.Math.BooleanVectors.firstTrueIndex({if i>firstnonzero[r] then nu[r,i] <> 0 else false for i in 1:nL}) for r in 1:nR};

  MassFraction[nX] Xrednonzero;

  Real[nL,nL] J;
  Real[nL] f;
  Real[nL] x1;
  Real[nL] x2;

  Real[nR,nL] B;
  Real[nR] b;

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

algorithm

  //ensure no species amount to be zero for numerical issues
  for i in 1:nX loop
    Xrednonzero[i] :=max(Modelica.Constants.eps, Xred[i]);
  end for;

  // find initial guess x1 for Newton solver
  for r in 1:nR loop
    B[r,secondnonzero[r]] :=1;
    b[r] :=0;
  end for;
  Xfull_in := ElectrolyteMedia.Media.Common.equalityLeastSquares(
    transpose(lambda_mass),
    Xred,
    B,
    b);
  // nfull_in with positive elements only to fulfill transpose(lambda_mass)*mfull_in = Xred;
  if sum(Xfull_in) > 0 then
    mfull_in :=Xfull_in;
    nfull_in :=mfull_in ./ MMX;
    for i in 1:nF loop
      x1[i] :=max(init, nfull_in[i]);
    end for;
  elseif sum(Xinit)> 0 then
  // ninit from initial guess provided as input (Xinit)
    scale :=sum(transpose(lambda_mass)*Xinit);
    minit :=Xinit*scale;
    Yinit :=Functions.calc_Y(Xinit);
    ninit :=minit./MMX;
    for i in 1:nF loop
      x1[i] :=max(init, ninit[i]);
    end for;
  else
  // else generic initial guess
    x1 :=fill(1/nF, nF);
  end if;

  // Newton algorithm
  solutionfound :=false;
  while not solutionfound loop
    f :=Initialization.calc_f.LE_Tp(x1,Xrednonzero,T,p);
    J :=Initialization.calc_J.LE_Tp(x1,T,p);
    x2 := Initialization.NewtonStep(x1,f,J,nL);
    if Modelica.Math.Vectors.norm(f,Modelica.Constants.inf) < eps then
      solutionfound :=true;
    end if;
    x1 :=x2;
  end while;

  Yfull :=x2/sum(x2);
  Xfull :=Functions.calc_X(Yfull);

end X_pTXred;
