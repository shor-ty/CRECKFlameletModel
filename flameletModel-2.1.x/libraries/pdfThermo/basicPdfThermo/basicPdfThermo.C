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

\*---------------------------------------------------------------------------*/

#include "basicPdfThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(basicPdfThermo, 0);
    defineRunTimeSelectionTable(basicPdfThermo, fvMesh);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicPdfThermo::basicPdfThermo(const fvMesh& mesh)
:
    basicThermo(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::basicPdfThermo::~basicPdfThermo()
{}


// Virtual stuff
Foam::volScalarField& Foam::basicPdfThermo::csi()
{
    notImplemented("basicPdfThermo::csi()");
    return const_cast<volScalarField&>(volScalarField::null());
}

const Foam::volScalarField& Foam::basicPdfThermo::csi() const
{
    notImplemented("basicPdfThermo::csi() const");
    return volScalarField::null();
}

Foam::volScalarField& Foam::basicPdfThermo::csiv2()
{
    notImplemented("basicPdfThermo::csiv2()");
    return const_cast<volScalarField&>(volScalarField::null());
}

const Foam::volScalarField& Foam::basicPdfThermo::csiv2() const
{
    notImplemented("basicPdfThermo::csiv2() const");
    return volScalarField::null();
}

Foam::volScalarField& Foam::basicPdfThermo::chi_st()
{
    notImplemented("basicPdfThermo::chi_st()");
    return const_cast<volScalarField&>(volScalarField::null());
}

const Foam::volScalarField& Foam::basicPdfThermo::chi_st() const
{
    notImplemented("basicPdfThermo::chi_st() const");
    return volScalarField::null();
}

Foam::volScalarField& Foam::basicPdfThermo::H()
{
    notImplemented("basicPdfThermo::H()");
    return const_cast<volScalarField&>(volScalarField::null());
}

const Foam::volScalarField& Foam::basicPdfThermo::H() const
{
    notImplemented("basicPdfThermo::H() const");
    return volScalarField::null();
}

Foam::volScalarField& Foam::basicPdfThermo::as()
{
    notImplemented("basicPdfThermo::as()");
    return const_cast<volScalarField&>(volScalarField::null());
}

const Foam::volScalarField& Foam::basicPdfThermo::as() const
{
    notImplemented("basicPdfThermo::as() const");
    return volScalarField::null();
}

// ************************************************************************* //
