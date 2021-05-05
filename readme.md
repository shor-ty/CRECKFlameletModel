##  Rebuild of the steady-state laminar Flamelet Model for OpenFOAM v8

The libraries you are using here are developed by Alberto Cuoci and his team (CRECK Modeling group). For more information have a look at the official website http://creckmodeling.chem.polimi.it/

## Supported OpenFOAM versions
 - OpenFOAM v8
 - OpenFOAM v7
 - OpenFOAM 4.x
 - OpenFOAM 2.3.x
 - OpenFOAM 2.2.x
 - OpenFOAM 2.1.x

## Introduction
Tobias worked with that flamelet-tool over 8 month during his Master thesis in 2012 and analysed the whole thermo model and the flamelet extraction libraries. For that Tobias made some validations and compared the results with ANSYS-CFX and measurements. There are some "bugs" in the original version which are removed in this version. Furthermore, you can use the SIMPLEC algorithm for solving steady-state combustions and the PIMPLE algorithm for transient calculations.

The library only contains a steady-state solver. If you create the transient one, feel free to make a push request.


## Discussion
A discussion can be followed here: http://www.cfd-online.com/Forums/openfoam-programming-development/99645-libopensmoke.html

## Compiling

In the documentation folder you will find a compilation instruction which should be valid for this version too. However, you should
first try the following steps given below. If you do have any problems, feel free to mail to community@Holzmann-cfd.com. The »master« branch represents the actual supported versions (here OpenFOAM v8).

```bash
cd $FOAM_RUN
cd ../
git clone https://github.com/shor-ty/CRECKFlameletModel.git CRECKFlameletModel
cd CRECKFlameletModel
```
Now, depending on your OpenFOAM version, you have to replace »master« with the branch that fits to you.
```bash
git checkout master
cd thermophysicalModels/flameletExtraction/common
wmake libso
cd ../flamelets/turbulent/
wmake libso
cd ../../../basic/
wmake libso
cd ../../applications/solvers/combustion/flameletSimpleFoam
wmake
cd ../../../../
cp -r tutorials $FOAM_RUN/flameletTutorials
```

## Changes | Features | Documentation
+ Have a look into the documentation folder

## Important | Validation

+ the rebuilt of the flamelet model for OpenFOAM v7 and v8 was not validated. Only the necessary class changes were performed in order to use the library for newer versions. Nevertheless, one should first countercheck the SANDIA flame tutorial with the measurment data to ensure correctness.

## Notice | Warranty
+ Not tested. No warranty of results and accuraty.

## Older versions
+ the modified version of 7.x can be found inside this repository
+ the modified version of 4.x can be found inside this repository
+ the modified version of 2.3.x can be found inside this repository
+ the modified version of 2.2.x can be found inside this repository
+ the modified version of 2.1.x can be found inside this repository
+ the origin version of 2.1.x can be found here: https://github.com/wyldckat/libOpenSMOKE
+ the origin version of 1.7.x can be found here: https://github.com/wyldckat/libOpenSMOKE


