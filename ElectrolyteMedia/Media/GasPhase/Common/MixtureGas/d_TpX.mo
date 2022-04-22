within ElectrolyteMedia.Media.GasPhase.Common.MixtureGas;
function d_TpX "Calculates density of gas"
  extends Modelica.Icons.Function;

  input SI.Temperature T "Temperature";
  input SI.AbsolutePressure p "Pressure";
  input MassFraction[:] X;

  output SI.Density d "Temperature";

protected
  MassFraction[nX] Xfull = if size(X,1) == nX then X else cat(1,X,{1-sum(X)});
    package Internal
    "Solve h(data,T) for T with given h (use only indirectly via temperature_phX)"
    extends ElectrolyteMedia.Media.Common.OneNonLinearEquation;

     redeclare function extends f_nonlinear
     algorithm
     y :=Functions.calc_p(
               T_,
               x,
               X_);
     end f_nonlinear;

    // Dummy definition has to be added for current Dymola
      redeclare function extends solve
      end solve;
    end Internal;

    MoleFraction[:] Y = massToMoleFractions(X,MMX);
    SI.MolarVolume v_min;
    SI.MolarVolume v_max;
    SI.Density d_max;
    SI.Density d_min;

    SI.Pressure p_temp;
    SI.Pressure p_min;
    SI.Density d_temp;

    Real i;
    parameter Real stepsize = 1e-3;

    Boolean otherzero = false;
    Boolean nootherzero = false;

algorithm
  if GasModel == Media.Common.Types.GasModel.Ideal then
    d :=p/(X*data[:].R*T);

  elseif GasModel == Media.Common.Types.GasModel.PengRobinson then
    v_min:=
      Functions.PengRobinson.Molar.calc_b(
       Y)*1.00001;
    v_max :=1e4;
    d_max:=Y*MMX/v_min;
    d_min:=Y*MMX/v_max;
    d := Internal.solve(
      p,
      d_min,
      d_max,
      0,
      T,
      0,
      X);
    // check if there is also solution for lower density. For that purpose, span grid and check if zero crossing within grid. If so, find zero of reduced interval
      i :=stepsize;

      p_min :=Functions.calc_p(
        T,
        d_min,
        X) - p;

      while not otherzero and not nootherzero loop
        d_temp :=d - (d - d_min)*i;
        p_temp :=Functions.calc_p(
          T,
          d_temp,
          X) - p;
        i :=i + stepsize;
        if p_temp > 0 and p_min < 0 or p_temp < 0 and p_min > 0 then
          otherzero :=true;
        end if;
        if i > 1 then
          nootherzero :=true;
        end if;
      end while;
     if otherzero then
     d :=Internal.solve(
       p,
       d_min,
       d_temp,
       0,
       T,
       0,
       X);
     end if;
  end if;

  annotation(experiment(Tolerance=1e-08));
end d_TpX;
