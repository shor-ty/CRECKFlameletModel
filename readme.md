##  Rebuild version of the flamelet model for OpenFOAM 4.x

The libraries you are using here are developed by Alberto Cuoci and his team. For more information have a look at the official website http://creckmodeling.chem.polimi.it/

## Introduction

Rebuild libraries for OpenFOAM-4.x. More information of introduction can be read in older version

## Discussion
A discussion can be followed here: http://www.cfd-online.com/Forums/openfoam-programming-development/99645-libopensmoke.html

## Compiling

Please do the following steps to compile the library and the application

```bash
cd $FOAM_RUN
cd ../
git clone https://github.com/shor-ty/flameletModel.git
cd flameletModel
git checkout flameletModel-2.3.x
cd flameletModel-4.x/thermophysicalModels/flameletExtraction/common
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

+ the rebuild flamelet model for OpenFOAM 2.3.x is not validated but the changes which had made in the thermodynamic should not change the results 
+ same library for as flameletModel-2.2.x


## Notice | Warranty
+ Not tested. No warranty of results and accuraty.

## Older versions
+ the modified version of 2.3.x can be found within this repository
+ the modified version of 2.2.x can be found within this repository
+ the modified version of 2.1.x can be found within this repository
+ the origin version of 2.1.x can be found here: https://github.com/wyldckat/libOpenSMOKE
+ the origin version of 1.7.x can be found here: https://github.com/wyldckat/libOpenSMOKE



