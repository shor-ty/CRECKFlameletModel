##  Rebuild of the steady-state laminar Flamelet Model for OpenFOAM 7.x

The libraries you are using here are developed by Alberto Cuoci and his team (CRECK Modeling group). For more information have a look at the official website http://creckmodeling.chem.polimi.it/

## Introduction

Tobias worked with that flamelet-tool over 8 month during his Master thesis in 2012 and analysed the whole thermo model and the flamelet extraction libraries. For that Tobias made some validations and compared the results with ANSYS-CFX and measurements. There are some "bugs" in the original version which are removed in this version. Furthermore, you can use the SIMPLEC algorithm for solving steady-state combustions and the PIMPLE algorithm for transient calculations.


## Discussion
A discussion can be followed here: http://www.cfd-online.com/Forums/openfoam-programming-development/99645-libopensmoke.html

## Compiling

IN the documentation folder you will find a compilation instruction which should be valid for this version too. However, you should
first try the following steps given above. If you do have any problems, feel free to mail to community@Holzmann-cfd.com.

```bash
cd $FOAM_RUN
cd ../
git clone https://github.com/shor-ty/CRECKFlameletModel.git CRECKFlameletModel
cd CRECKFlameletModel
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

+ the rebuild flamelet model for OpenFOAM 7.x is not validated but the changes which had made in the thermodynamic should not change the results
+ same library for as flameletModel-4.x

## Notice | Warranty
+ Not tested. No warranty of results and accuraty.

## Older versions
+ the modified version of 4.x can be found within this repository
+ the modified version of 2.3.x can be found within this repository
+ the modified version of 2.2.x can be found within this repository
+ the modified version of 2.1.x can be found within this repository
+ the origin version of 2.1.x can be found here: https://github.com/wyldckat/libOpenSMOKE
+ the origin version of 1.7.x can be found here: https://github.com/wyldckat/libOpenSMOKE


