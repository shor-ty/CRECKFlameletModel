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
    Foam::hPdfThermo

Description
    Enthalpy for a mixture based on compressibility

SourceFiles
    hPdfThermo.C

\*---------------------------------------------------------------------------*/

#ifndef hPdfThermo_H
#define hPdfThermo_H

#include "basicPdfThermo.H"
#include "basicMixture.H"

#include "OpenSMOKE_PDF_NonAdiabaticFlamelet_Library.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class hPdfThermo Declaration
\*---------------------------------------------------------------------------*/

template<class MixtureType>
class hPdfThermo
:
    public basicPdfThermo,
    public MixtureType
{
private:


	volScalarField h_;

	volScalarField csi_;
	volScalarField csiv2_;
	volScalarField H_;

	volScalarField density_reynolds_;
	volScalarField chi_st_;
	volScalarField phiH_;
	volScalarField as_;		
	volScalarField mu_favre_;	
	volScalarField alpha_favre_;
	PtrList<volScalarField> omega_;
				
	OpenSMOKE_PDF_NonAdiabaticFlamelet_Library 	flamelets_library;

	double HFuel; 
	double HOxidizer;
	label counter;
	label propertyUpdate;
	label counter_mass_fractions;
	label massFractionsUpdate;
	label nDirac;
	Switch adiabaticMode;

	std::vector<int> patch_type_T;
	std::vector<int> patch_type_H;
	std::vector<int> patch_type_csi;

private:

	hPdfThermo(const hPdfThermo<MixtureType>&);

	void calculate();
	void update();
	void updateMassFractions();
	void ErrorMessage(const string message);


public:

        //- Read thermophysicalProperties dictionary
        virtual bool read();
		
	//- Mixture fraction
	//  Non-const access allowed for transport equations
	virtual volScalarField& csi()
	{
		return csi_;
	}

	//- Mixture fraction
	virtual const volScalarField& csi() const
	{
		return csi_;
	}

	//- Variance of Mixture fraction
	//  Non-const access allowed for transport equations
	virtual volScalarField& csiv2()
	{
		return csiv2_;
	}

	//- Variance of Mixture fraction
	virtual const volScalarField& csiv2() const
	{
		return csiv2_;
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

public:

    //- Runtime type information
    TypeName("hPdfThermo");


    // Constructors

        //- Construct from mesh
        hPdfThermo(const fvMesh&);


    //- Destructor
    virtual ~hPdfThermo();


    // Member functions

        //- Return the compostion of the mixture
        virtual basicMixture& composition()
        {
            return *this;
        }

        //- Return the compostion of the mixture
        virtual const basicMixture& composition() const
        {
            return *this;
        }

        //- Update properties
        virtual void correct();


        // Access to thermodynamic state variables

            //- Enthalpy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& h()
            {
                return h_;
            }

            //- Enthalpy [J/kg]
            virtual const volScalarField& h() const
            {
                return h_;
            }


        // Fields derived from thermodynamic state variables

            //- Enthalpy for cell-set [J/kg]
            virtual tmp<scalarField> h
            (
                const scalarField& T,
                const labelList& cells
            ) const;

            //- Enthalpy for patch [J/kg]
            virtual tmp<scalarField> h
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant pressure for patch [J/kg/K]
            virtual tmp<scalarField> Cp
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant pressure [J/kg/K]
            virtual tmp<volScalarField> Cp() const;

            //- Heat capacity at constant volume for patch [J/kg/K]
            virtual tmp<scalarField> Cv
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant volume [J/kg/K]
            virtual tmp<volScalarField> Cv() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
#   include "hPdfThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
