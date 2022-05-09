# `ElectrolyteMedia` in Modelica

The `ElectrolyteMedia` library is based on the Modelica Standard Library `Media` package and allows the simulation of customized aqueous electrolyte systems under consideration of chemical equilibrium. The library comprises multiple packages for the simulation of gas, liquid, and solid phases and combinations thereof. The object-oriented package extends the standardized Modelica structure for calculating intensive thermodynamic properties of a fluid medium incorporating thermodynamic models. The equilibrium conditions are then incorporated in the `BaseProperties` model of the considered system. The intensive thermodynamic properties calculated in the `BaseProperties` models find use in component models comprising mass and energy balances. Modelica provides generic models containing mass and energy balances in the `Fluid` library of the Modelica Standard Library. Alternatively, user-specific models with customized mass and energy balances can be implemented and interfaced with the `BaseProperties` models to calculate intensive thermodynamic properties and the equilibrium composition.

## How to use the `ElectrolyteMedia` library

The package contains a [user manual](ElectrolyteMedia/Media/MediumUsage.mo) that provides the following:

- Overview of implemented `Media` packages with considered thermodynamic models
- Introduction of reaction invariants used as independent mass fractions and how they are used within the `Media` packages
- Guide how to set up a user-specific medium by using the user interface and providing parameters of thermodynamic models.

The framework provides multiple [test models](ElectrolyteMedia/Tests/) that show the capabilities of the `Media` packages. We refer to a [titration model](ElectrolyteMedia/Tests/LiquidPhase/NaOH_HCL_titration.mo) that simulates the titration of an NaOH solution with a HCl solution. The model extends from a [generic titration model](ElectrolyteMedia/Tests/LiquidPhase/partial_Titration.mo). Both models incorporate detailed explanations of the modeling steps for reproduciton and/or modeling user-specific aqueous electrolyte systems.

We will shortly add a reference to a submitted manuscript that discusses the underlying methods for the simulation of aqueous electroylte systems under consideration of chemical equilibrium.


