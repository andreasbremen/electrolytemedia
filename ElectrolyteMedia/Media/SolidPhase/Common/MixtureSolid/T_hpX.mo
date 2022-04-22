within ElectrolyteMedia.Media.SolidPhase.Common.MixtureSolid;
function T_hpX "Compute temperature from specific enthalpy"
  extends Modelica.Icons.Function;
  input SpecificEnthalpy h "Specific enthalpy";
  input AbsolutePressure p;
  input MassFraction[nX] X;
  output Temperature T "Temperature";

protected
    package Internal
    "Solve h(data,T) for T with given h (use only indirectly via temperature_phX)"
    extends Media.Common.OneNonLinearEquation;

      redeclare function extends f_nonlinear
      algorithm
      y := Functions.calc_h(x,p_,X_);//h_TpX(x,p_,X_);
      end f_nonlinear;

    // Dummy definition has to be added for current Dymola
      redeclare function extends solve
      end solve;
    end Internal;

algorithm
  T := Internal.solve(h, 273.15, 573.15, p,1,1, X);
end T_hpX;
