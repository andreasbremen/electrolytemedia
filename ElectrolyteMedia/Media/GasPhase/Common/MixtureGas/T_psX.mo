within ElectrolyteMedia.Media.GasPhase.Common.MixtureGas;
function T_psX "Calculates temperature of gas"
  extends Modelica.Icons.Function;

  input SI.AbsolutePressure p "Pressure";
  input SI.SpecificEntropy s "Specific Enthalpy";
  input MassFraction[:] X;

  output SI.Temperature T "Temperature";

protected
  MassFraction[nX] Xfull = if size(X,1) == nX then X else cat(1,X,{1-sum(X)});
    package Internal
    "Solve h(data,T) for T with given h (use only indirectly via temperature_phX)"
    extends ElectrolyteMedia.Media.Common.OneNonLinearEquation;

     redeclare function extends f_nonlinear
     algorithm
     y :=s_TpX(x,p_,X_);
     end f_nonlinear;

    // Dummy definition has to be added for current Dymola
      redeclare function extends solve
      end solve;
    end Internal;

algorithm
  T := Internal.solve(
    s,
    273.15,
    573.15,
    p,
    0,
    0,
    X);

  annotation(smoothOrder(normallyConstant = data)=20);
end T_psX;
