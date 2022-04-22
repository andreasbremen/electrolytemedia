within ElectrolyteMedia.Media.GasLiquidPhase.Common.SaturatedLiquid;
function T_Kwp
  "Calculates boiling temperature according to equilibrium constant of water and pressure"
  extends Modelica.Icons.Function;

  input Real Kw "Equilibrium constant";
  input SI.AbsolutePressure p "Pressure";

  output SI.Temperature T "Temperature";

//   MassFraction[nG] Xfull = if size(X,1) == nG then X else cat(1,X,{1-sum(X)});
protected
    package Internal
    "Solve h(data,T) for T with given h (use only indirectly via temperature_phX)"
    extends Media.Common.OneNonLinearEquation;

     redeclare function extends f_nonlinear
     algorithm
     y :=calc_Kw(x,p_);
     //Functions.GasFunctions.calc_p(T_,x,X_);
     end f_nonlinear;

    // Dummy definition has to be added for current Dymola
      redeclare function extends solve
      end solve;
    end Internal;

protected
  parameter SI.Temperature T_min = 260;
  parameter SI.Temperature T_max = 600;

algorithm
  T :=Internal.solve(
    Kw,
    T_min,
    T_max,
    p,
    0,
    0,
    fill(1/(nG+nL),nG+nL));
end T_Kwp;
