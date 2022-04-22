within ElectrolyteMedia.Tests.GasLiquidPhase.MixtureLiquid;
model Test_volume
  Modelica.Fluid.Vessels.ClosedVolume volume(
    redeclare package Medium =
        Media.GasLiquidPhase.MixtureLiquids.ExampleMedium,
    energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial,
    massDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial,
    p_start=200000,
    T_start=298.15,
    X_start=volume.Medium.reference_X,
    use_portsData=false,
    V=1,
    medium(
      preferredMediumStates=false,
      Tstart=volume.T_start,
      pstart=volume.p_start,
      Xredstart=volume.X_start),
    nPorts=1) annotation (Placement(transformation(extent={{-10,50},{10,70}})));

  Modelica.Fluid.Sources.MassFlowSource_T boundary(
    redeclare package Medium =
        Media.GasLiquidPhase.MixtureLiquids.ExampleMedium,
    use_m_flow_in=true,
    T=323.15,
    X=volume.Medium.reference_X,
    nPorts=1,
    medium(
      Tstart=volume.T_start,
      pstart=volume.p_start,
      Xredstart=boundary.X))
    annotation (Placement(transformation(extent={{-50,30},{-30,50}})));
  inner Modelica.Fluid.System system
    annotation (Placement(transformation(extent={{-40,60},{-20,80}})));

  Modelica.Blocks.Sources.Constant mass_flow(k=10)
    annotation (Placement(transformation(extent={{-90,38},{-70,58}})));
equation

  connect(volume.ports[1], boundary.ports[1])
    annotation (Line(points={{0,50},{0,40},{-30,40}},  color={0,127,255}));
  connect(boundary.m_flow_in, mass_flow.y)
    annotation (Line(points={{-50,48},{-69,48}}, color={0,0,127}));
 annotation (Line(points={{-85,4},
          {-72,4},{-72,-12},{-60,-12}}, color={0,0,127}),
              Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)),
    experiment(
      StopTime=100,
      __Dymola_NumberOfIntervals=100,
      Tolerance=1e-08,
      __Dymola_Algorithm="Lsodar"));
end Test_volume;
