within ElectrolyteMedia.Media.GasLiquidPhase.Common.MixtureLiquid;
function X_pTinert
  "Gas-liquid and liquid dissociation equilibrium with output mole fractions of both liquid and gas phase separately. Fixed mole fractions of inert gas species and fixed molalities of inert liquid species"
   input SI.Pressure p;
   input SI.Temperature T;
   input Real[nG+nL] inert;
   input MoleFraction[nG+nL] xinit = zeros(nF);
   output MoleFraction[nG+nL] x_ "Mass fraction for each gas and liquid phase";
protected
  Real[nG + nL,nG + nL] J;
  Real[nG + nL] f;
  Real[2,nG + nL] x = ones(2,nG + nL);
  Real eps = 1e-4;
  Boolean solutionfound = false;
algorithm
  if sum(xinit) > 0 then
    x[1,:] :=xinit;
  else
    for i in 1:nG loop
      x[1,i] :=1/nG;
    end for;
    x[1,1+nG:nG+nL-1] :=fill(1e-7,nL-1);
    x[1,nG+nL] :=1;
  end if;
  while not solutionfound loop
    f :=Initialization.calc_f.GLE_Tp(
      x[1, :],
      inert,
      T,
      p);
    J :=Initialization.calc_J.GLE_Tp(
      x[1, :],
      inert,
      T,
      p);
    x[2, :] :=Initialization.NewtonStep(
      x[1, :],
      f,
      J,
      nG + nL);

    if Modelica.Math.Vectors.norm((x[2,:]-x[1,:])./(x[2,:]),Modelica.Constants.inf)< eps then
      solutionfound :=true;
    end if;
    x[1,:] :=x[2,:];
  end while;
  x_ :=x[1, :];
end X_pTinert;
