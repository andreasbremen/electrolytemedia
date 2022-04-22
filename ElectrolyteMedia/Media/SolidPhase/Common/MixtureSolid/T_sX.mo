within ElectrolyteMedia.Media.SolidPhase.Common.MixtureSolid;
function T_sX "Compute temperature from pressure and specific entropy"
  extends Modelica.Icons.Function;
  input SpecificEntropy s "Specific entropy";
  input MassFraction[nX] X;
  output Temperature T "Temperature";

protected
    package Internal
    "Solve h(data,T) for T with given h (use only indirectly via temperature_phX)"
    extends Media.Common.OneNonLinearEquation;

      redeclare function extends f_nonlinear
      algorithm
      y := Functions.calc_s(x,X_);
      end f_nonlinear;

    // Dummy definition has to be added for current Dymola
      redeclare function extends solve
      end solve;
    end Internal;

algorithm
  T := Internal.solve(s,273.15,573.15,1.0e5,1,1,X);
end T_sX;
