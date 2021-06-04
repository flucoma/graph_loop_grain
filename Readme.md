# Graph-based audio looping and granulation

This repository contains the source code for three experimental objects for audio playback, based on the self-similarity graph derived from time-frequency analysis. Links between spectral frames of high similarity can be seen as 'wormholes', that can be used to create alternative playback paths. The graph is also used to find clusters of similar frames, detect onsets and find common repetitive periods.
Three objects are included for Max and SuperCollider: FluidGraphlLoop (Max: fluid.graphloop~) implements a looper that finds good looping points based on similarity. FluidGraphGrain (Max: fluid.graphgrain~) implements a granular synthesis algorithm that concatenates spectral frames within a cluster. FluidGraphPlay (Max: fluid.graphplay~) implements a stochastic playback algorithm that allows jumping at random locations based on similarity and a minimum segment length.

The objects are based on the [Fluid Corpus Manipulation Library](https://flucoma.org)

#  Building for Max

## Dependencies

See the [Max repository of the Fluid Corpus Manipulation Library](https://github.com/flucoma/flucoma-max) for the dependencies. In particular, you need the [Max SDK](https://github.com/Cycling74/max-sdk) (>= 7.3.3). All other dependencies are downloaded automatically during the build process.

## Compilation

Starting from the `source/max` folder of this repository:
```
mkdir build && cd build
cmake -DMAX_SDK_PATH=<location_of_your_Max_SDK> ..
make install
```

This should result in a `dist` folder with the Max objects and demo packages.


#  Building for SuperCollider

## Dependencies

See the [SuperCollider repository of the Fluid Corpus Manipulation Library](https://github.com/flucoma/flucoma-sc) for the dependencies. In particular, you need the  [SuperCollider Source Code](https://github.com/supercollider/supercollider). All other dependencies are downloaded automatically during the build process.

## Compilation

Starting from the `source/supercollider` folder of this repository:
```
mkdir build && cd build
cmake -DSC_PATH=<location_of_your_SC_source> ..
make install
```

This should result in a `dist` folder with the plugins, classes and help files.
