within ElectrolyteMedia.Media.GasLiquidPhase.Common.SaturatedLiquid;
function p_TdX "calculates pressure of liquid mixture"
  extends Modelica.Icons.Function;

  input SI.Temperature T "Temperature";
  input SI.Density d "Pressure";
  input MassFraction[:] X;

  output SI.Pressure p;

protected
    package Internal
    "Solve h(data,T) for T with given h (use only indirectly via temperature_phX)"
    extends Media.Common.OneNonLinearEquation;

     redeclare function extends f_nonlinear
     algorithm
      y :=Functions.calc_d(
             T_,
             x,
             X_);
     end f_nonlinear;

    // Dummy definition has to be added for current Dymola
      redeclare function extends solve
      end solve;
    end Internal;

protected
    SI.Pressure pmin = Modelica.Media.Water.IF97_Utilities.BaseIF97.Basic.psat(T)*0.3;

     SI.SpecificVolume v = 1/d;

algorithm
  p := Internal.solve(
    d,
    pmin,
    100e6,
    0,
    T,
    0,
    X);   //pmin,

  annotation(smoothOrder(normallyConstant = data)=20, experiment(Tolerance=1e-08));
end p_TdX;
