/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::pdfFlameletThermo

Description
    Energy for a mixture based on density

SourceFiles
    pdfFlameletThermo.C

\*---------------------------------------------------------------------------*/

#ifndef pdfFlameletThermo_H
#define pdfFlameletThermo_H

#include "flameletThermo.H"
#include "heThermo.H"
#include "OpenSMOKE_PDF_NonAdiabaticFlamelet_Library.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class pdfFlameletThermo Declaration
\*---------------------------------------------------------------------------*/

template<class BasicPsiThermo, class MixtureType>
class pdfFlameletThermo
:
    public heThermo<BasicPsiThermo, MixtureType>
{
    // Private Member Functions

        //- Calculate the thermo variables
        void calculate();

        //- Update hte thermo variables by extracting from LUT
        void update();

        //- Update the species mass fraction by extracting from LUT
        void updateMassFractions();

        //- Error message
        void errorMessage(const string message);

        //- Info message
        void infoMessage() const;

        //- Construct as copy (not implemented)
        pdfFlameletThermo
        (
            const pdfFlameletThermo<BasicPsiThermo,
            MixtureType>&
        );


public:

    //- Runtime type information
    TypeName("pdfFlameletThermo");


    // Constructors

        //- Construct from mesh and phase name
        pdfFlameletThermo
        (
            const fvMesh&,
            const word& phaseName
        );


    //- Destructor
    virtual ~pdfFlameletThermo();


    // Member functions

        //- Update properties
        virtual void correct();

private:

    // Variables

        //- Mixture fraction
        volScalarField Z_;

        //- Variance of mixture fraction
        volScalarField Zvar_;

        //- Scalar dissipation rate
        volScalarField chi_st_;

        //- Enthalpy of mixture 
        volScalarField H_;

        //- Enthalpy defect due to non-adiabatic system
        volScalarField defect_;

        //- PLank absorbtions coeff (for radiation)
        volScalarField as_;

        //- Reynolds density
        volScalarField density_reynolds_;

        //- Dynamic viscosity (favre)
        volScalarField mu_favre_;

        //- Thermal diffusivity
        volScalarField alpha_favre_;
        

        //- Species mass fraction
        PtrList<volScalarField> omega_;

        //- Adiabatic enthalpy of fuel
        scalar HFuel;

        //- Adiabatic enthalpy of oxidizer
        scalar HOxidizer;

        //- Variable for set thermo variable extraction
        label counter;

        //- variable for set mass fraction extraction
        label counter_mass_fractions;

        //- User defined numbers for update
        label propertyUpdate;
        label massFractionsUpdate;

        //label nDirac;
        Switch adiabaticMode;
        Switch showFlamelet;
        Switch showFlameletLibrary;

        //- stuff (modifyed by Holzmann)
        List<int> patch_type_T;
        List<int> patch_type_H;
        List<int> patch_type_Z;


private:

    //- Flamelet class
    OpenSMOKE_PDF_NonAdiabaticFlamelet_Library flamelets_library;


public:

    //- Mixture fraction

        //  Non-const access allowed for transport equations
        virtual volScalarField& Z()
        {
            return Z_;
        }

        //- Mixture fraction
        virtual const volScalarField& Z() const
        {
            return Z_;
        }

        //- Variance of Mixture fraction
        //  Non-const access allowed for transport equations
        virtual volScalarField& Zvar()
        {
            return Zvar_;
        }

        //- Variance of Mixture fraction
        virtual const volScalarField& Zvar() const
        {
            return Zvar_;
        }

        //- Scalar dissipation rate [1/s]
        //  Non-const access allowed for transport equations
        virtual volScalarField& chi_st()
        {
            return chi_st_;
        }

        //- Scalar dissipation rate [1/s]
        virtual const volScalarField& chi_st() const
        {
            return chi_st_;
        }

        //- Enthalpy [J/kg]
        //  Non-const access allowed for transport equations
        virtual volScalarField& H()
        {
            return H_;
        }

        //- Enthalpy [J/kg]
        virtual const volScalarField& H() const
        {
            return H_;
        }

        //- Absorption coefficient [1/m]
        //  Non-const access allowed for transport equations
        virtual volScalarField& as()
        {
            return as_;
        }

        //- Absorption coefficient [1/m]
        virtual const volScalarField& as() const
        {
            return as_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "pdfFlameletThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
