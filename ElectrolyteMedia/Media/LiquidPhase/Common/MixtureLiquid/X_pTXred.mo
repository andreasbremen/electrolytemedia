within ElectrolyteMedia.Media.LiquidPhase.Common.MixtureLiquid;
function X_pTXred
  "Liquid dissociation equilibrium calculated via Newton method with homotopy continuation method"
  input SI.Pressure p;
  input SI.Temperature T;
  input MassFraction[nX] Xred;
  input MassFraction[nF] Xinit = zeros(nF);
  output Real[nF] x_;

protected
  parameter Real init = 1e-10;
  parameter Real eps = 1e-6;

  MassFraction[nX] Xrednonzero;
  MassFraction[nX] Xred_;

  Real[nF,nF] J;
  Real[nF] f;
  Real[nF] x1;
  Real[nF] x2;

  SI.MoleFraction[nF] Yfull;
  SI.MassFraction[nF] Xfull;

  SI.MoleFraction[nF] Yinit;
  SI.AmountOfSubstance[nF] ninit;
  SI.Mass[nF] minit;
  Real scale;

  Boolean solutionfound = false;

  Real[nX] Xs;
algorithm
  //ensure no species amount to be zero for numerical issues
  for i in 1:nX loop
    Xrednonzero[i] :=max(Modelica.Constants.eps, Xred[i]);
  end for;

  //convert reduced mass fraction basis
  Xs :=Modelica.Math.Matrices.solve(transpose(lambda_mass_id[nR+1:nF, :]), Xrednonzero);
  Xfull :=cat(1,zeros(nF - nX),Xs);
  Xfull :=P_to_orig*Xfull;

  Xred_ :=transpose(lambda_mass_id_orig)*Xfull;

  // find initial guess x1 for Newton solver
  if sum(Xinit)> 0 then
  // ninit from initial guess provided as input (Xinit)
     scale :=sum(transpose(lambda_mass_id_orig)*Xinit);
     minit :=Xinit*scale;
     ninit :=minit./MMX;
     for i in 1:nF loop
       x1[i] :=max(init, ninit[i]);
     end for;
  else
  // else generic initial guess
    for i in 1:nF loop
      x1[i] :=max(init, Xfull[i]/MMX[i]);
    end for;
  end if;

  // Newton algorithm
  solutionfound :=false;
  x1 :=P_to_id*x1;
  while not solutionfound loop
    f :=Initialization.calc_f.LE_Tp(x1,Xred_,T,p);
    J :=Initialization.calc_J.LE_Tp(x1,T,p);
    x2 := Initialization.NewtonStep(x1,f,J,nF);

    // check convergence
    if Modelica.Math.Vectors.norm(f,Modelica.Constants.inf) < eps then
      solutionfound:=true;
    end if;
    x1 :=x2;
  end while;
  x2 :=P_to_orig*x2;
  Yfull :=x2[1:nF]/sum(x2[1:nF]);
  Xfull :=moleToMassFractions(Yfull, MMX);

  x_ :=Xfull;

end X_pTXred;
