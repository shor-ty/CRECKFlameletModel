/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "flameletThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(flameletThermo, 0);
    defineRunTimeSelectionTable(flameletThermo, fvMesh);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::flameletThermo::flameletThermo(const fvMesh& mesh, const word& phaseName)
:
    fluidThermo(mesh, phaseName),

    psi_
    (
        IOobject
        (
            phasePropertyName("thermo:psi"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(0, -2, 2, 0, 0)
    ),

    mu_
    (
        IOobject
        (
            phasePropertyName("thermo:mu"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(1, -1, -1, 0, 0)
    )
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::flameletThermo> Foam::flameletThermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    return basicThermo::New<flameletThermo>(mesh, phaseName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::flameletThermo::~flameletThermo()
{}


// * * * * * * * * * * * * * * * Virtual Stuff  * * * * * * * * * * * * * * * //

Foam::volScalarField& Foam::flameletThermo::Z()
{
    notImplemented("flameletThermo::Z()");
    return const_cast<volScalarField&>(volScalarField::null());
}

const Foam::volScalarField& Foam::flameletThermo::Z() const
{
    notImplemented("flameletThermo::Z() const");
    return volScalarField::null();
}

Foam::volScalarField& Foam::flameletThermo::Zvar()
{
    notImplemented("flameletThermo::Zvar()");
    return const_cast<volScalarField&>(volScalarField::null());
}

const Foam::volScalarField& Foam::flameletThermo::Zvar() const
{
    notImplemented("flameletThermo::Zvar() const");
    return volScalarField::null();
}

Foam::volScalarField& Foam::flameletThermo::chi_st()
{
    notImplemented("flameletThermo::chi_st()");
    return const_cast<volScalarField&>(volScalarField::null());
}

const Foam::volScalarField& Foam::flameletThermo::chi_st() const
{
    notImplemented("flameletThermo::chi_st() const");
    return volScalarField::null();
}

Foam::volScalarField& Foam::flameletThermo::H()
{
    notImplemented("flameletThermo::H()");
    return const_cast<volScalarField&>(volScalarField::null());
}

const Foam::volScalarField& Foam::flameletThermo::H() const
{
    notImplemented("flameletThermo::H() const");
    return volScalarField::null();
}

Foam::volScalarField& Foam::flameletThermo::as()
{
    notImplemented("flameletThermo::as()");
    return const_cast<volScalarField&>(volScalarField::null());
}

const Foam::volScalarField& Foam::flameletThermo::as() const
{
    notImplemented("flameletThermo::as() const");
    return volScalarField::null();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::flameletThermo::rho() const
{
    return p_*psi();
}

const Foam::volScalarField& Foam::flameletThermo::psi() const
{
    return psi_;
}


const Foam::volScalarField& Foam::flameletThermo::mu() const
{
    return mu_;
}


// ************************************************************************* //
