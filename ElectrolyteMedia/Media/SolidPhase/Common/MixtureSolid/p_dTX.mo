within ElectrolyteMedia.Media.SolidPhase.Common.MixtureSolid;
function p_dTX "Calculates pressure from density and temperature"

  extends Modelica.Icons.Function;
  input Temperature T "Temperature";
  input Density d;
  input MassFraction[nX] X;
  output AbsolutePressure p "Specific enthalpy";

protected
    package Internal
    "Solve h(data,T) for T with given h (use only indirectly via temperature_phX)"
    extends Media.Common.OneNonLinearEquation;

      redeclare function extends f_nonlinear
      algorithm
      y :=Functions.calc_d(T_,x,X_);
      end f_nonlinear;

    // Dummy definition has to be added for current Dymola
      redeclare function extends solve
      end solve;
    end Internal;

algorithm
  p := Internal.solve(d, 1e5, 1000e5, 1e5,T,d, X);

end p_dTX;
