within ElectrolyteMedia.Media.SolidPhase.MixtureSolid;
package ExampleMedium
  extends Common.MixtureSolid(userInterface(
    nX=1,
    substanceNames={"NaCl"},
    data=Common.SolidData.halite,
    mediumName="halite",
    Tstart=298.15,
    pstart=100000,
    useSolidMassFraction=true));

end ExampleMedium;
