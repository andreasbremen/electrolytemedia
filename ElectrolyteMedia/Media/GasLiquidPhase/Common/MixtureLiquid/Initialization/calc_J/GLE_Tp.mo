within ElectrolyteMedia.Media.GasLiquidPhase.Common.MixtureLiquid.Initialization.calc_J;
function GLE_Tp
  "Gas-liquid and dissociation equilibrium with constraints on inert solute molalities and moles of gas phase"
  input Real[nG+nL] x "mole fractions of gas and liquid phase for each phase separately";
  input Real[nG+nL] inert "mole fraction of inert gas phase species and molalities of inert liquid phase species";
  input SI.Temperature T;
  input SI.Pressure p;
  output Real[nG + nL,nG + nL] J;

protected
  Boolean[nG+nL] isinert = {sum(abs(nu[:,i])) < Modelica.Constants.eps for i in 1:nG+nL};
  Integer index;
algorithm
  for i in  1:nG + nL loop
    //isopotential
    for r in 1:nR loop
      if i < nG+1 then
        J[r,i] := if nu[r,i]<> 0 then nu[r, i]/x[i] else 0;
      elseif i == nG+nL then
         J[r,i] :=-sum({nu[r, i_] for i_ in 1 + nG:nG + nL - 1})/x[i] + nu[r, i]/x[i];
      else J[r,i] :=if abs(nu[r,i]) > 0 then nu[r, i]/x[i] else 0;
      end if;
    end for;

    //nR+1: closure condition gas phase
    //nR+2: closure condition liquid phase
    //nR+3: charge balance liquid phase
    if i < nG+1 then
      J[nR+1,i] :=1;
      J[nR+2,i] :=0;
      J[nR+3,i] :=0;
    elseif i == nG+nL then
      J[nR+1,i] :=0;
      J[nR+2,i] :=1;
      J[nR+3,i] :=0;
    else
      J[nR+1,i] :=0;
      J[nR+2,i] :=1;
      J[nR+3,i] :=datal[i-nG].z;
    end if;
  end for;

  //constraints on inert solute molalities and mole fraction of gas phase
  index :=4;
  for i in 1:nG+nL-1 loop
    if isinert[i] then
      if i < nG+1 then
        J[nR+index,i] :=1000;
      else
        J[nR+index,i] :=1;
        J[nR+index,nG+nL] :=-MH2O*inert[i];
      end if;
      index :=index + 1;
       if index > nG+nL-nR then
         break;
       end if;
    end if;
  end for;
end GLE_Tp;
