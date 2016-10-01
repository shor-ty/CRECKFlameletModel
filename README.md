# README #

This README  includes the steps which are necessary to get this application up and running.

### OpenFOAM Versions ###
* Generated for version 
** 2.1.x **
** 2.2.x **
** 2.3.x **
** 4.x **

### What is this repository for? ###
* This model contains the thermodynamic flamelet model and a binary file for building Look-Up-Tables for the model
* combustion model for RANS simulation for non-premixed turbulent flow fields (especially for free stream flames)
* It is a combustion model without solving chemistry during fluid dynamic simulation
* Chemistry calculation is done as a pre-processing and stored in Look-Up-Tables which are used for the simulation
* Developed by [CRECK Modelling](www.creckmodeling.chem.polimi.it)
* Rebuild and modified by [Holzmann-cfd](www.holzmann-cfd.de)

### Features | Changes | Documentation###
* In the repository folder you will find PDF's with all necessary stuff and validations

### Validation ###
* You find some validation hints in the folder too
* Additionally you can check out some validation at [cfd-online.com](http://www.cfd-online.com/Forums/openfoam-verification-validation/115150-validation-flamelet-model-rebuild-2-2-x.html)

### Discussion ###
* There is a Thread on cfd-online where you can find interesting things: [libOpenSMOKE](http://www.cfd-online.com/Forums/openfoam-programming-development/99645-libopensmoke.html)

### Older Versions ###
* There is also an old version for 1.7.x and 2.1.x (not modified libraries) available: [github - wyldkcat](https://github.com/wyldckat/libOpenSMOKE)

### How do I get set up? ###
* Feel free to compile it where ever you want, but normally its nice to have a fixed folder for _user compiled stuff_
* Make a new folder
> mkdir -p $FOAM_RUN/../OpenFOAM_extensions
* Switch to the new folder
> cd $FOAM_RUN/../OpenFOAM_extensions
* Clone the repository to the new folder
> git clone git@bitbucket.org:shor-ty/flameletmodel.git
* Switch to the repository directory
> cd flameletModel
* Switch to the branch you need

* for version 2.1.x use the following command
> git checkout flameletModel-2.1.x
* for version 2.2.x use the following command
> git checkout flameletModel-2.2.x
* for version 2.3.x use the following command
> git checkout flameletModel-2.3.x
* for version 4.x use the following command
> git checkout flameletModel-4.x
* After that pull it
> git pull
* Now you have your libraries available and you can start compiling
* Start with the common libraries (x.x.x » stand for your version)
> cd flameletModel-x.x.x/thermophysicalModels/flameletExtraction/common
>> wmake libso
* Go on with the flamelet library
> cd ../flamelets/turbulent/
>> wmake libso
* Hence you finished that build the thermodynamic library
> cd ../../../basic/
>> wmake libso
* Now you can build the flameletSimpleFoam solver
> cd ../../applications/solvers/combustion/flameletSimpleFoam
>> wmake
* At least go back
> cd ../../../../
* And copy the tutorials to your run folder
> cp -r tutorials $FOAM_RUN/flameletTutorials
* Finished

### Contribution guidelines ###
* If you have questions, hints or any suggestions please email me to Tobias.Holzmann@Holzmann-cfd.de

### Other stuff ###
* Thank to Alexander Vakrushev for introduce me to bitbucket and git program, also for OpenFOAM in general
* Thanks to Bruno Santos for rebuild the model for my master thesis, and also a introduction to git and branches
