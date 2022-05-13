within ElectrolyteMedia.Tests;
package SolidLiquidPhase

  model Test_volume_Modelica
    Modelica.Fluid.Vessels.ClosedVolume volume(
      redeclare package Medium =
          Media.SolidLiquidPhase.MixtureLiquids.ExampleMedium,
      energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial,
      massDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial,
      p_start=1000000,
      T_start=323.15,
      X_start=volume.Medium.reference_X,
      use_portsData=false,
      V=1,
      medium(
        preferredMediumStates=false,
        Tstart=volume.T_start,
        pstart=volume.p_start,
        Xredstart=volume.X_start),
      nPorts=1)
      annotation (Placement(transformation(extent={{-10,0},{10,20}})));

    Modelica.Fluid.Sources.MassFlowSource_T boundary(
      redeclare package Medium =
          Media.SolidLiquidPhase.MixtureLiquids.ExampleMedium,
      use_m_flow_in=false,
      m_flow=1,
      T=323.15,
      X=volume.Medium.reference_X,
      nPorts=1,
      medium(
        Tstart=volume.T_start,
        pstart=volume.p_start,
        Xredstart=boundary.X))
      annotation (Placement(transformation(extent={{-60,-30},{-40,-10}})));
    inner Modelica.Fluid.System system
      annotation (Placement(transformation(extent={{60,60},{80,80}})));

  equation

    connect(volume.ports[1], boundary.ports[1])
      annotation (Line(points={{0,0},{0,-20},{-40,-20}}, color={0,127,255}));
   annotation (Line(points={{-85,4},
            {-72,4},{-72,-12},{-60,-12}}, color={0,0,127}),
                Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)),
      experiment(
        __Dymola_NumberOfIntervals=100,
        Tolerance=1e-06,
        __Dymola_Algorithm="Lsodar"));
  end Test_volume_Modelica;

  partial model Base
    "Test method to generate outputs from given inputs"

    //Declare Medium
    replaceable package Medium =
        Media.SolidLiquidPhase.Common.MixtureLiquid;

    //Declare ThermodynamicState
    //To be specified
  //    Real[Medium.nL-1] m_in;//mofified
    parameter Real[Medium.nL-1] m0;
  //    SI.MassFraction[Medium.ns] Xs_in;//modified
    parameter SI.MassFraction[Medium.ns] Xs0;
  //    Real liquidFraction; //modified
    parameter Real liquidFraction0;
    //Given
    parameter Temperature T0=300;
    parameter AbsolutePressure p0=1e5;
  //   SI.MassFraction[Medium.nL] Xl_in = Medium.Functions.LiquidFunctions.calc_Xfromm(m_in);
    parameter SI.MassFraction[Medium.nL] Xl0 = Medium.Functions.LiquidFunctions.calc_Xfromm(m0);
  //   SI.MassFraction[Medium.nF] X_in = cat(1,(1-liquidFraction)*Xs_in,liquidFraction*Xl_in);//modified
    parameter SI.MassFraction[Medium.nF] X0 = cat(1,(1-liquidFraction0)*Xs0,liquidFraction0*Xl0);
  //   SI.MassFraction[Medium.nX] Xred_in = transpose(Medium.lambda_mass)*X_in/sum(transpose(Medium.lambda_mass)*X_in);//modified
    parameter SI.MassFraction[Medium.nX] Xred0 = transpose(Medium.lambda_mass)*X0/sum(transpose(Medium.lambda_mass)*X0);

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)),
      experiment(Tolerance=1e-08));
  end Base;

  model setState_pTX
    extends Base(
      redeclare package Medium =
          Media.SolidLiquidPhase.MixtureLiquids.NaCl_dis,
      m0={0.01,0,0.01,0},
      Xs0={1},
      liquidFraction0=0.95);

    Medium.ThermodynamicState state;
    Medium.BaseProperties medium(Xredstart = Xred0,Tstart=T0,pstart=p0);

    SpecificEntropy s_calc;
    SpecificEnthalpy h_calc;
    Density d_calc;

  equation

    s_calc = Medium.specificEntropy(state);
    h_calc = Medium.specificEnthalpy(state);
    d_calc = Medium.density(state);

    state = Medium.setState_pTX(p0, T0, Xred0);

    medium.T = state.T;
    medium.p = state.p;
    medium.X[1:Medium.nX-1] = state.X[1:Medium.nX-1];

    annotation (experiment(Tolerance=1e-06));
  end setState_pTX;

  model setState_phX
    extends Base(
      redeclare package Medium =
          Media.SolidLiquidPhase.MixtureLiquids.NaCl_dis,
      m0={0.01,0,0.01,0},
      Xs0={1},
      liquidFraction0=0.95);

    parameter SI.SpecificEnthalpy h0 = -15412027;

    Medium.ThermodynamicState state;
    Medium.BaseProperties medium(Xredstart = Xred0,Tstart=T0,pstart=p0);

  equation

    state = Medium.setState_phX(p0, h0, Xred0);

    medium.T = state.T;
    medium.p = state.p;
    medium.X[1:Medium.nX-1] = state.X[1:Medium.nX-1];

    annotation (experiment(Tolerance=1e-06));
  end setState_phX;

  model setState_psX
    extends Base(
      redeclare package Medium =
          Media.SolidLiquidPhase.MixtureLiquids.NaCl_dis,
      m0={0.01,0,0.01,0},
      Xs0={1},
      liquidFraction0=0.95);

    parameter SI.SpecificEntropy s0 = 3809.0745;

    Medium.ThermodynamicState state;
    Medium.BaseProperties medium(Xredstart = Xred0,pstart=p0,Tstart=T0);

  equation

    state = Medium.setState_psX(p0, s0, Xred0);

    medium.T = state.T;
    medium.p = state.p;
    medium.X[1:Medium.nX-1] = state.X[1:Medium.nX-1];

    annotation (experiment(Tolerance=1e-06));
  end setState_psX;

  model setState_dTX
    extends Base(
      redeclare package Medium =
          Media.SolidLiquidPhase.MixtureLiquids.NaCl_dis,
      m0={0.01,0,0.01,0},
      Xs0={1},
      liquidFraction0=0.95);

    parameter SI.Density d0 = 1034.352;

    Medium.ThermodynamicState state;
    Medium.BaseProperties medium(Xredstart = Xred0,Tstart=T0,pstart=p0);

  equation

    state = Medium.setState_dTX(d0, T0, Xred0);

    medium.T = state.T;
    medium.p = state.p;
    medium.X[1:Medium.nX-1] = state.X[1:Medium.nX-1];

    annotation (experiment(Tolerance=1e-06));
  end setState_dTX;

  model ReactiveTransportModel

    replaceable package Medium =
        Media.SolidLiquidPhase.Common.MixtureLiquid        annotation(choicesAllMatching=true);

    parameter SI.Temperature T = 333.15;
    parameter SI.Pressure p = 100e5;

    parameter Real[Medium.nX] Xred = Medium.userInterface.refX;

    parameter Integer N = 21 "Number of finite elements";
    parameter Integer nX = Medium.nX;
    parameter Integer ns = Medium.ns;
    parameter Integer nL = Medium.nL;
    parameter Real[:,:] lambda_mass = Medium.lambda_mass;

    parameter SI.Length L = 0.2;
    parameter SI.Length deltaL = L/(N-1);
    parameter SI.Length B = 1;
    parameter SI.Length H = 1;
    parameter SI.Volume V = B*H*deltaL;

    parameter SI.VolumeFraction Phil0 = 0.5;
    parameter SI.Velocity u = 1.2e-5;

    parameter Real[nL-1] ml_in = fill(1e-9,nL-1);
    parameter SI.MassFraction[nL] Xl_in = Medium.Functions.LiquidFunctions.calc_Xfromm(ml_in);
    parameter SI.MassFraction[ns] Xs_in = fill(1/ns,ns);
    parameter SI.MassFraction lfrac_in = 1;
    parameter SI.MassFraction[ns+nL] Xfull_in = cat(1,(1-lfrac_in)*Xs_in,lfrac_in*Xl_in);
    parameter SI.Density d_in = Medium.Functions.calc_d(T,p,Xfull_in);
    parameter SI.Density[ns+nL] d_in_i = d_in*Xfull_in;
    parameter SI.Density dl_in = Medium.Functions.LiquidFunctions.calc_d(T,p,Xl_in);

    parameter SI.MassFraction[ns + nL] Xfullstart=Medium.X_pTXred(p,T,Medium.reference_X);

    parameter SI.Density d_start = Medium.Functions.calc_d(T,p,medium_n[1].Xfullstart);
    parameter SI.Mass m_start = d_start*V;

    Medium.BaseProperties[N] medium_n(each Tstart=T,each pstart=p, each Xfullstart = Xfullstart);
    SI.Mass[nX,N] mred_n(stateSelect = StateSelect.always);
    SI.MassFraction[ns+nL,N] Xfull_flow;
    SI.Velocity[N] v_n_flow(each start=u/Phil0);
    SI.Density[ns+nL,N] d_i;
    SI.MolarDensity[ns+nL,N] dm_i;
    SI.Mass[N] mtot_n;
    SI.Mass[ns+nL,N] m_n;
    SI.MassFlowRate[ns+nL,N] m_in;
    SI.MassFlowRate[ns+nL,N] m_out;
  initial equation
    for n in 1:N loop
       mred_n[:,n] = transpose(lambda_mass)*Medium.userInterface.refXfull*d_start*B*H*deltaL;
    end for;

  equation

    // constitutive equations
    for n in 1:N loop
      Xfull_flow[:,n] = cat(1,zeros(ns),medium_n[n].Xl);
      for i in 1:Medium.nF loop
        d_i[i,n] = medium_n[n].d*medium_n[n].Xfull[i];
        dm_i[i,n] = d_i[i,n]/Medium.MMX[i];
      end for;
      // volume constraint
      V = mtot_n[n]/medium_n[n].d;
      m_n[:,n] = mtot_n[n]*medium_n[n].Xfull;

      // flows in and out of each volume element
      m_out[:,n] = v_n_flow[n]*B*H*medium_n[n].Phil*medium_n[n].dl*Xfull_flow[:,n];
      if n > 1 then
        m_in[:,n] = m_out[:,n-1];
      else
        m_in[:,n] = u*B*H*dl_in*Xfull_in;
      end if;
    end for;

    // differential euqations
    for n in 1:N loop
        for j in 1:nX loop
          der(mred_n[j,n]) = sum(lambda_mass[i,j]*(m_in[i,n]-m_out[i,n]) for i in 1:ns+nL);
        end for;
      der(mtot_n[n]) = sum(m_in[:,n]-m_out[:,n]);
    end for;

    // connect to medium
    for n in 1:N loop
      medium_n[n].T = T;
      medium_n[n].p = p;
      medium_n[n].X[1:nX-1] = mred_n[1:nX-1,n]/sum(mred_n[:,n]);
    end for;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)),
      experiment(
        StopTime=60000,
        Tolerance=1e-06,
        __Dymola_Algorithm="Lsodar"));
  end ReactiveTransportModel;

  model ReactiveTransportSimulation

    ReactiveTransportModel reactiveTransport(
      redeclare package Medium =
          Media.SolidLiquidPhase.MixtureLiquids.ReactiveTransport,
      T=333.15,
      p=100e5,
      ml_in={0.75,0,0,0.05,0.01,0.9,1.02,0,0,0},
      lfrac_in=1)
      annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
    annotation (experiment(
      StopTime=60000,
      Tolerance=1e-06,
      __Dymola_Algorithm="Lsodar"));
  end ReactiveTransportSimulation;
end SolidLiquidPhase;
