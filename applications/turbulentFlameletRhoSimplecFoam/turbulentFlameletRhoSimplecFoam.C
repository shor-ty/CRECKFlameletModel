/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Application
    rhoSimplecFoam

Description
    Steady-state SIMPLEC solver for laminar or turbulent RANS flow of
    compressible fluids.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "basicPdfThermo.H"
#include "hPdfThermo.H"
#include "RASModel.H"
#include "mixedFvPatchFields.H"
#include "bound.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readGravitationalAcceleration.H"
    #include "readMassFlowProperties.H"

    simpleControl simple(mesh);

    #include "createFields.H"
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Velocity-pressure-enthalpy SIMPLEC corrector
        {
            #include "UEqn.H"
            #include "pEqn.H"
            #include "hEqn.H"
	    #include "csiEqn.H"
        }

        turbulence->correct();

        runTime.write();

	Info<< "- - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
	Info<< "rho   min|max: \t" << min(rho).value() << "\t" << max(rho).value() << endl;
	Info<< "csi   min|max: \t" << min(csi).value() << "\t\t" << max(csi).value() << endl;
	Info<< "csiv  min|max: \t" << min(csiv2).value() << "\t\t" << max(csiv2).value() << endl;
	Info<< "H     min|max: \t" << min(H).value() << "\t" << max(H).value() << endl;
	Info<< "- - - - - - - - - - - - - - - - - - - - - - - - -" << endl;

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
