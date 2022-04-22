within ElectrolyteMedia.Media.GasLiquidPhase.Common.MixtureLiquid;
function X_pTXred
  "Gas-liquid and liquid dissociation equilibrium with T, p and reduced mass fractions (sum(Xred)=1) as input. Output is quasi-homogeneous mass fraction of gas and liquid species."
  input SI.Pressure p;
  input SI.Temperature T;
  input MassFraction[nX] Xred;
  input MassFraction[nF] Xinit = zeros(nF);
  output MassFraction[nF] x_;

protected
  parameter Real eps = 1e-5;

  Real[nF,nF] J;
  Real[nF] f;
  Real[nF] x1;
  Real[nF] x2;

  Boolean solutionfound = false;

  parameter Integer[nR] firstnonzero = {Modelica.Math.BooleanVectors.firstTrueIndex({nu[r,i] <> 0 for i in 1:nF}) for r in 1:nR};

  Boolean[nF] isinert = {sum(abs(nu[:,i])) < Modelica.Constants.eps for i in 1:nG+nL};
  Real lfrac;
  Real[nF] inert;

  Real[nR,nF] B;
  Real[nR] b;

  SI.MassFraction[nG] Xg;
  SI.MoleFraction[nG] Yg;
  SI.MassFraction[nL] Xl;
  SI.MoleFraction[nL] Yl;
  Real[nL-1] ml;

  SI.MassFraction[nF] Xfull_in;

  SI.AmountOfSubstance[nF] Xnotinert;
  SI.AmountOfSubstance[nX] Xrednotinert;

  Real[nX] Xred_;

  SI.MoleFraction[nF] Yfull;
  SI.MassFraction[nF] Xfull;
  SI.MoleFraction[nF] Yinit;

  Integer index;
algorithm

  // find initial guess nfull_in with positive elements only and 1kg in total
  for r in 1:nR loop
    B[r,firstnonzero[r]] :=1;
    b[r] :=0;
  end for;
  Xfull_in :=Modelica.Math.Matrices.equalityLeastSquares(transpose(lambda_mass), Xred,B,b);
  Xfull_in :=Xfull_in/sum(Xfull_in);

  // initial guess for Newton solver
  if sum(Xinit)> 0 then
    Yinit :=Functions.calc_Y(Xinit);
    Yinit[1:nG] :=Yinit[1:nG]/sum(Yinit[1:nG]);
    Yinit[1+nG:nF] :=Yinit[1 + nG:nF]/sum(Yinit[1 + nG:nF]);
  end if;

  for i in 1:nG+nL loop
    if not isinert[i] then
      Xnotinert[i] :=Xfull_in[i];
    end if;
  end for;

  //find index for calculation of liquid fraction with non-inert species
  Xrednotinert :=transpose(lambda_mass)*Xnotinert;
  if index < Modelica.Constants.eps then
    index :=Modelica.Math.BooleanVectors.firstTrueIndex({Xrednotinert[i] > 0 for i in 1:nX});
  end if;
  Xred_ :=transpose(lambda_mass)*Xfull_in;

  // specify inert: fixed mole fractions of inert gas species and fixed molalities of inert liquid species
  Xg :=Xfull_in[1:nG]/(sum(Xfull_in[1:nG]) + Modelica.Constants.eps);
  Yg :=if sum(Xg) > 0 then Functions.GasFunctions.calc_Y(Xg) else zeros(nG);
  Xl :=Xfull_in[1+nG:nF]/sum(Xfull_in[1+nG:nF]);
  ml :=Functions.LiquidFunctions.calc_mfromX(Xl);
  for i in 1:nG+nL loop
    if isinert[i] then
       if i < nG+1 then
         inert[i] :=min(0.9,Yg[i]);// leads to errors if Yg[i] = 1
       else
        inert[i] :=ml[i-nG];
       end if;
    end if;
  end for;

  while not solutionfound loop
    // calculate mole fractions of gas phase and liquid phase
    Yfull :=X_pTinert(p,T,inert,Yinit);

    // determine liquid fraction lfrac from noninert species
    x1[1:nG] :=Functions.GasFunctions.calc_X(Yfull[1:nG]);
    x1[1+nG:nG+nL] :=Functions.LiquidFunctions.calc_X(Yfull[1+nG:nF]);
    lfrac :=(Xred_[index] - x1[1:nG]*lambda_mass[1:nG, index])/(x1[1 + nG:nG + nL]*
      lambda_mass[1 + nG:nG + nL, index] - x1[1:nG]*lambda_mass[1:nG, index]);
    for i in 1:nG+nL loop
      if i < nG+1 then
        x1[i] :=x1[i]*(1 - lfrac);
      else
        x1[i] :=x1[ i]*lfrac;
      end if;
    end for;

    // lfrac needs to be between 0 and 1
    assert(lfrac < 1, "Input mass fraction is not in two phase region --> only liquid phase!",AssertionLevel.error);
    assert(lfrac > 0, "Input mass fraction is not in two phase region --> only gas phase!",AssertionLevel.error);

    // define break criterion
    if Modelica.Math.Vectors.norm((x2-x1)./(x2),Modelica.Constants.inf)< eps then
      solutionfound :=true;
    end if;

    //store iteration result in x2
    x2 :=x1;

    // update inert species for next iteration (mole fractions for gaseous species and molalities for liquid species)
    Xg :=x1[1:nG];
    for i in 1:nG loop
      if isinert[i] then
        Xg[i] :=Xfull_in[i];
      end if;
    end for;
    Xg :=Xg/sum(Xg);
    Yg :=Functions.GasFunctions.calc_Y(Xg);

    Xl :=x1[1 + nG:nF];
    Yl :=Functions.LiquidFunctions.calc_Y(Xl);
    for i in 1:nL loop
      if isinert[i+nG] then
        Xl[i] :=Xfull_in[i + nG];
      end if;
    end for;
    Xl :=Xl/sum(Xl);
    ml :=Functions.LiquidFunctions.calc_mfromX(Xl);

    // fixed mole fractions of inert gas species and fixed molalities of inert liquid species
    for i in 1:nG+nL-1 loop
      if isinert[i] then
         if i < nG+1 then
           inert[i] :=Yg[i];
         else
           inert[i] :=ml[i-nG];
         end if;
      end if;
    end for;

    //provide start value for next iteration
    Yinit :=cat(1,Yg,Yl);
  end while;

  x_ :=x1;

end X_pTXred;
