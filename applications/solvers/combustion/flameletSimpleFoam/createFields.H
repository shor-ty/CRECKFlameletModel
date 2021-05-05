Info<< "Reading thermophysical properties\n" << endl;

autoPtr<flameletThermo> pThermo
(
    flameletThermo::New(mesh)
);
flameletThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

Info<< "Preparing scalar fields from flameletThermo\n"
    << "  + thermo.p()\n"
    << "  + thermo.Z()\n"
    << "  + thermo.Zvar()\n"
    << "  + thermo.chi_st()\n"
    << "  + thermo.H()\n"
    << "  + thermo.as()\n"
    << endl;

volScalarField& p = thermo.p();
volScalarField& Z = thermo.Z();
volScalarField& Zvar = thermo.Zvar();
volScalarField& chi_st = thermo.chi_st();
volScalarField& H = thermo.H();
volScalarField& as = thermo.as();

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "compressibleCreatePhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, simple.dict(), pRefCell, pRefValue);

mesh.setFluxRequired(p.name());

dimensionedScalar rhoMax
(
    dimensionedScalar::lookupOrDefault
    (
        "rhoMax",
        simple.dict(),
        dimDensity,
        GREAT
    )
);

dimensionedScalar rhoMin
(
    dimensionedScalar::lookupOrDefault
    (
        "rhoMin",
        simple.dict(),
        dimDensity,
        0
    )
);

Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::momentumTransportModel> turbulence
(
    compressible::momentumTransportModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

Info<< "Reading flameletProperties dictionary\n" << endl;
IOdictionary flameletsProperties
(
    IOobject
    (
        "flameletsProperties",
        U.time().constant(),
        U.db(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

dimensionedScalar initialMass = fvc::domainIntegrate(rho);

#include "createMRF.H"

//- Monitoring
Switch monitoring(flameletsProperties.lookup("monitoring"));

//- Semi-empirical constants in mixture fraction variance equation
scalar Cg = readScalar(flameletsProperties.lookup("Cg"));
scalar Cd = readScalar(flameletsProperties.lookup("Cd"));

//- Semi-empirical constant for scalar dissipation rate evaluation
scalar Cx = readScalar(flameletsProperties.lookup("Cx"));

//- Prandtl number for mixture fraction, variance of mixture fraction
//  and enthalpy diffusion terms
dimensionedScalar sigmat
    (
        "sigmat",
        dimensionSet(0,0,0,0,0,0,0),
        readScalar(flameletsProperties.lookup("sigmat"))
    );

//- Minimum values (to be used for scalar dissipation rate evaluation)
dimensionedScalar k_small(flameletsProperties.lookup("k_small"));
dimensionedScalar epsilon_small
    (
        flameletsProperties.lookup("epsilon_small")
    );

//- Environment temperature (non adiabatic consitions)
dimensionedScalar Tenv (flameletsProperties.lookup("Tenv"));

//- Variance of mixture fraction:
//      on = transport equation 
//      off = algebraic equation
Switch varianceEquation(flameletsProperties.lookup("varianceEquation"));

//- Adiabatic simulation:
//      on = adiabatic
//      off = non adiabatic
Switch adiabaticMode(flameletsProperties.lookup("adiabaticMode"));
Switch radiationMode(flameletsProperties.lookup("radiationMode"));

bool switchQoff = false;

Info<< "Preparing field Qrad (radiative heat transfer)\n" << endl;
volScalarField Qrad
(
    IOobject
    (
        "Qrad",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("qradiation", dimensionSet(1,-1,-3,0,0,0,0) , 0.0)
);

dimensionedScalar sigmaSB
    (
        "stefan_bolzmann",
        dimensionSet(1,0,-3,-4,0,0,0),
        5.6704e-8
    );


