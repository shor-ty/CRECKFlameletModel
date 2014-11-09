# README #

This README  includes the steps which are necessary to get this application up and running.

### OpenFOAM Versions ###
* Generated for Version 2.3.x

### What is this repository for? ###
* This post-processing tool calculates the total kinetic energy in our domain for any time folder.
* It can be used for compressible and incompressible fluids
* Version 1.0
* Developed by [Holzmann-cfd](https://holzmann-cfd.de)

### How do I get set up? ###
* Feel free to compile it where ever you want, but normally its nice to have a fixed folder for _user compiled stuff_
* Make a new folder
> mkdir -p $FOAM_RUN/../OpenFOAM_extensions
* Switch to the new folder
> cd $FOAM_RUN/../OpenFOAM_extensions
* Clone the repository to the new folder
> git clone https://shor-ty@bitbucket.org/shor-ty/totalkineticenergy.git
* Switch to the repository directory
> cd totalKineticEnergy
* Compile the application
> wmake
* Finished

### Contribution guidelines ###
* If you have questions, hints or any suggestions please email me to Tobias.Holzmann@Holzmann-cfd.de

### Other stuff ###
* Thank to Alexander Vakrushev for a short introduction to the C++ OpenFOAM functions
