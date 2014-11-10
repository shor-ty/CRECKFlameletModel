/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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
    Foam::basicPdfThermo

Description
    Basic thermodynamic properties based on compressibility

SourceFiles
    basicPdfThermo.C
    newBasicPdfThermo.C

\*---------------------------------------------------------------------------*/

#ifndef basicPdfThermo_H
#define basicPdfThermo_H

#include "volFields.H"
#include "basicThermo.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class basicPdfThermo Declaration
\*---------------------------------------------------------------------------*/

class basicPdfThermo
:
    public basicThermo
{

protected:

    // Protected member functions

        //- Construct as copy (not implemented)
        basicPdfThermo(const basicPdfThermo&);


public:

    //- Runtime type information
    TypeName("basicPdfThermo");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        basicPdfThermo,
        fvMesh,
        (const fvMesh& mesh),
        (mesh)
    );

    // Constructors

        //- Construct from mesh
        basicPdfThermo(const fvMesh&);


    //- Selector
    static autoPtr<basicPdfThermo> New(const fvMesh&);

    //- Destructor
    virtual ~basicPdfThermo();


    // Member functions

        // Virtual stuff
        virtual volScalarField& csi();
        virtual const volScalarField& csi() const;

        virtual volScalarField& csiv2();
        virtual const volScalarField& csiv2() const;

        virtual volScalarField& chi_st();
        virtual const volScalarField& chi_st() const;

        virtual volScalarField& H();
        virtual const volScalarField& H() const;

        virtual volScalarField& as();
        virtual const volScalarField& as() const;
    
        // Fields derived from thermodynamic state variables

            //- Density [kg/m^3] - uses current value of pressure
            virtual tmp<volScalarField> rho() const
            {
                return p_*psi();
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
