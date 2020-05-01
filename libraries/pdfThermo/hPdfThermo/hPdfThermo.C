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

#include "hPdfThermo.H"
#include "fixedValueFvPatchFields.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class MixtureType>
Foam::hPdfThermo<MixtureType>::hPdfThermo(const fvMesh& mesh)
:
    basicPdfThermo(mesh),
    MixtureType(*this, mesh),

    h_
    (
        IOobject
        (
            "h",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(0, 2, -2, 0, 0),
        hBoundaryTypes()
    ), 
    
    csi_
    (
        IOobject
        (
            "csi",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    csiv2_
    (
        IOobject
        (
            "csiv2",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    
    H_
    (
        IOobject
        (
            "H",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),   

    density_reynolds_
    (
        IOobject
        (
	   "density_reynolds",
	   mesh.time().timeName(),
	   mesh,
	   IOobject::NO_READ,
	   IOobject::NO_WRITE
        ),
        mesh,
	dimensionedScalar("density",dimensionSet(1,-3,0,0,0,0,0) , 0.0)
    ),      
    
    chi_st_
    (
        IOobject
        (
            "chi_st",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("chi_st",dimensionSet(0,0,-1,0,0,0,0) , 0.0)
    ),  
    
    phiH_
    (
        IOobject
        (
            "phiH",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("phiH",dimensionSet(0,2,-2,0,0,0,0) , 0.0)
    ),
    
    as_
    (
        IOobject
        (
            "as",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("as",dimensionSet(0,-1,0,0,0,0,0) , 0.0)
   ),
   
    mu_favre_
    (
        IOobject
        (
            "mu_lam",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("mu_lam",dimensionSet(1,-1,-1,0,0,0,0) , 0.0)
   ),
   
    alpha_favre_
    (
        IOobject
        (
            "alpha_lam",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("alpha_lam",dimensionSet(0,2,-1,0,0,0,0) , 0.0)
   ),
    
   adiabaticMode(false)
{
	// Recognize fixed temperature boundaries
    	forAll(T_.boundaryField(), patchi)
    	{
		if (isA<fixedValueFvPatchScalarField>(T_.boundaryField()[patchi]))
		{
			Info << "Patch fixed Temperature " << mesh.boundary()[patchi].name() << endl;
			patch_type_T.push_back(1);
		}
		else
			patch_type_T.push_back(0);
    	}


	// Added for realize the inlets
	// Recognize fixed enthalpie boundaries 
    	forAll(H_.boundaryField(), patchi)
    	{
		if (isA<fixedValueFvPatchScalarField>(H_.boundaryField()[patchi]))
		{
			Info << "Patch fixed Enthalpy " << mesh.boundary()[patchi].name() << endl;
			patch_type_H.push_back(1);
		}
		else
		{
			patch_type_H.push_back(0);			
		}
    	}
	// Recognize inlets 
    	forAll(csi_.boundaryField(), patchi)
    	{
		if (isA<fixedValueFvPatchScalarField>(csi_.boundaryField()[patchi]))
		{
			Info << "Patch fixed mixture fraction (inlet " << mesh.boundary()[patchi].name() << endl;
			patch_type_csi.push_back(1);
		}
		else
		{
			patch_type_csi.push_back(0);			
		}
    	}

	// Dictionary
	IOdictionary flameletsProperties_
	(
		IOobject
		(
		    "flameletsProperties",
		    csi_.time().constant(),
		    csi_.db(),
		    IOobject::MUST_READ,
		    IOobject::NO_WRITE
		)
	);
    
    	string libraryPath = flameletsProperties_.lookup("libraryPath");
    	string chiPDF = flameletsProperties_.lookup("pdf");
    	string list_of_species = flameletsProperties_.lookup("species");

	Switch adiabaticMode_(flameletsProperties_.lookup("adiabaticMode"));
    	adiabaticMode 		= adiabaticMode_;

   	propertyUpdate 		= readLabel(flameletsProperties_.lookup("propertyUpdate"));
   	massFractionsUpdate 	= readLabel(flameletsProperties_.lookup("massFractionsUpdate"));

	// Starting counters
	counter = propertyUpdate;
	counter_mass_fractions = massFractionsUpdate;

	// Options
	flamelets_library.SetLibraryPath(libraryPath);
	flamelets_library.SetSpeciesToExtract(list_of_species);	
	
	// Initialize mass fraction fields
	omega_.setSize(flamelets_library.number_of_species()+1);
        for (int j=0;j<flamelets_library.number_of_species()+1;j++)
        {
		if(j < flamelets_library.number_of_species())
		{
                std::string name_of_species = "omega_" + flamelets_library.species()[j+1];
		omega_.set(     j, new volScalarField(  IOobject(name_of_species,
                                                                                        mesh.time().timeName(),
                                                                                        mesh,
                                                                                        IOobject::NO_READ,
                                                                                        IOobject::AUTO_WRITE),
                                                                                        mesh,
                                                                                        dimensionedScalar("zero", dimensionSet(0,0,0,0,0,0,0), 0.0)
                                                                                 )
                                        );
		}
        }
		
	// Set adiabatic mode
	if (adiabaticMode == true)
	{
		Info << "Adiabatic mode..." << endl;
		flamelets_library.SetAdiabaticMode();
	}
		
	// Scalar dissipation rate distribution
	if (chiPDF == "logNormal")
	{
		Info << "Log-Normal Distribution for scalar dissipation rate..." << endl;
		
		scalar chi_sigma(readScalar(flameletsProperties_.lookup("sigma")));
    		label  chi_points(readScalar(flameletsProperties_.lookup("points")));
    
		flamelets_library.UnsetExcludeColdFlamelets();
		flamelets_library.SetLogNormalChiDistribution(chi_sigma, chi_points);
	}
	else if (chiPDF == "dirac")
	{
		Info << "Delta-Dirac Distribution for scalar dissipation rate..." << endl;
	}
	
	// Read flamelets
	flamelets_library.Read();
	flamelets_library.Summary();
	

	// Update basic fields
	Info << "Initialize basic fields... " << endl;
	HOxidizer = flamelets_library.enthalpy_f_oxidizer();
    	HFuel = flamelets_library.enthalpy_f_fuel();
	update();

	// Initialize enthalpy field
	Info << "Initialize enthalpy field (field)... " << endl;
    	scalarField& hCells = h_.internalField();
    	const scalarField& HCells = H_.internalField();
    	forAll(hCells, celli)
    	{
		hCells[celli] = HCells[celli];
    	}
    	forAll(h_.boundaryField(), patchi)
    	{
		h_.boundaryField()[patchi] == h(T_.boundaryField()[patchi], patchi);
    	}
    	hBoundaryCorrection(h_);

	Info << "Calculate (first time)... " << endl;
	calculate();
    	psi_.oldTime();   // Switch on saving old time
    
    	Info << " Fuel  enthalpy:    " << HFuel 	<< " J/kg" << endl;
    	Info << " Oxid. enthalpy:    " << HOxidizer 	<< " J/kg" << endl;
    	Info << " Fuel  temperature: " << flamelets_library.temperature_f_fuel() << " K" << endl;
    	Info << " Oxid. temperature: " << flamelets_library.temperature_f_oxidizer() << " K" << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class MixtureType>
Foam::hPdfThermo<MixtureType>::~hPdfThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class MixtureType>
void Foam::hPdfThermo<MixtureType>::calculate()
{
	const scalarField& pCells	= p_.internalField();
	const scalarField& HCells 	= H_.internalField();
	const scalarField& csiCells 	= csi_.internalField();
	const scalarField& RhoReynolds 	= density_reynolds_.internalField();
	const scalarField& muFavre 	= mu_favre_.internalField();
	const scalarField& alphaFavre 	= alpha_favre_.internalField();

	scalarField& psiCells 		= psi_.internalField();
	scalarField& muCells 		= mu_.internalField();
	scalarField& alphaCells 	= alpha_.internalField();
	scalarField& phiHCells 		= phiH_.internalField();
	
	forAll(csiCells, celli)
	{
		psiCells[celli] 	= RhoReynolds[celli]/pCells[celli];
		muCells[celli] 		= muFavre[celli];
		alphaCells[celli]	= alphaFavre[celli];
		phiHCells[celli]	= HCells[celli] - (HOxidizer+csiCells[celli]*(HFuel-HOxidizer));
	}
    
	// Boundaries
	forAll(T_.boundaryField(), patchi)
	{
		const fvPatchScalarField& pp		= p_.boundaryField()[patchi];
		const fvPatchScalarField& pH 		= H_.boundaryField()[patchi];
		const fvPatchScalarField& pcsi 		= csi_.boundaryField()[patchi];
		const fvPatchScalarField& pRhoReynolds 	= density_reynolds_.boundaryField()[patchi];
		const fvPatchScalarField& pmuFavre 	= mu_favre_.boundaryField()[patchi];
		const fvPatchScalarField& palphaFavre 	= alpha_favre_.boundaryField()[patchi];

		fvPatchScalarField& pT 		= T_.boundaryField()[patchi];
		fvPatchScalarField& ppsi 	= psi_.boundaryField()[patchi];   
		fvPatchScalarField& pphiH 	= phiH_.boundaryField()[patchi];
		fvPatchScalarField& pmu 	= mu_.boundaryField()[patchi];
		fvPatchScalarField& palpha 	= alpha_.boundaryField()[patchi];


		if (pT.fixesValue())
		{
			forAll(pT, facei)
			{
				ppsi[facei] 	= pRhoReynolds[facei]/pp[facei]; 
				pmu[facei] 	= pmuFavre[facei];
				palpha[facei] 	= palphaFavre[facei];
				pphiH[facei]    = pH[facei] - (HOxidizer+pcsi[facei]*(HFuel-HOxidizer));
			}
		}
		else
		{
			forAll(pT, facei)
			{

				ppsi[facei] 	= pRhoReynolds[facei]/pp[facei]; 
				pmu[facei] 	= pmuFavre[facei];
				palpha[facei] 	= palphaFavre[facei];
				pphiH[facei]    = pH[facei] - (HOxidizer+pcsi[facei]*(HFuel-HOxidizer));
			}
		}
	}
}

template<class MixtureType>
void Foam::hPdfThermo<MixtureType>::correct()
{
	// force the saving of the old-time values
	psi_.oldTime();

	if (counter == propertyUpdate)
	{
		Info << "Updating look-up table extractions..." << endl;
		update();
		counter = 0;
	}

	if (counter_mass_fractions == massFractionsUpdate)
	{
		Info << "Updating mass fraction extractions..." << endl;
		updateMassFractions();
		counter_mass_fractions = 0;
	}

	calculate();

	counter++;
	counter_mass_fractions++;
}


template<class MixtureType>
Foam::tmp<Foam::scalarField> Foam::hPdfThermo<MixtureType>::h( const scalarField& T, const labelList& cells ) const
{
	Info << "pdfThermo h(const scalarField& T, const labelList& cells)...";

	tmp<scalarField> th(new scalarField(T.size()));
	scalarField& h = th();

	const scalarField& HCells = H_.internalField();
	forAll(T, celli)
	{
		//h[celli] = this->cellMixture(cells[celli]).H(T[celli]);
		h[celli] = HCells[celli];
	}
	
	Info << " end.." << endl;
	return th;
}


template<class MixtureType>
Foam::tmp<Foam::scalarField> Foam::hPdfThermo<MixtureType>::h( const scalarField& T, const label patchi ) const
{
	Info << "pdfThermo h(const scalarField& T, const label patchi)...";

	tmp<scalarField> th(new scalarField(T.size()));
	scalarField& h = th();

	 const fvPatchScalarField& pH = H_.boundaryField()[patchi];
	forAll(T, facei)
	{
		//h[facei] = this->patchFaceMixture(patchi, facei).H(T[facei]);
		h[facei] = pH[facei];
	}

	Info << " end.." << endl;
	return th;
}


template<class MixtureType>
Foam::tmp<Foam::scalarField> Foam::hPdfThermo<MixtureType>::Cp( const scalarField& T,const label patchi ) const
{
	Info << "pdfThermo Cp(const scalarField& T, const label patchi)...";

	tmp<scalarField> tCp(new scalarField(T.size()));
	scalarField& cp = tCp();

	forAll(T, facei)
	{
		cp[facei] = this->patchFaceMixture(patchi, facei).Cp(T[facei]);
	}

	Info << " end.." << endl;
	return tCp;
}


template<class MixtureType>
Foam::tmp<Foam::volScalarField> Foam::hPdfThermo<MixtureType>::Cp() const
{
	Info << "pdfThermo Cp()...";

	const fvMesh& mesh = T_.mesh();

	tmp<volScalarField> tCp
	(
		new volScalarField
		(
		    IOobject
		    (
			"Cp",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		    ),
		    mesh,
		    dimensionSet(0, 2, -2, -1, 0)
		)
	);

	volScalarField& cp = tCp();

	forAll(T_, celli)
	{
		cp[celli] = this->cellMixture(celli).Cp(T_[celli]);
	}

	forAll(T_.boundaryField(), patchi)
	{
		cp.boundaryField()[patchi] = Cp(T_.boundaryField()[patchi], patchi);
	}

	Info << " end.." << endl;
	return tCp;
}


template<class MixtureType>
Foam::tmp<Foam::scalarField> Foam::hPdfThermo<MixtureType>::Cv( const scalarField& T, const label patchi ) const
{
	Info << "pdfThermo Cv(const scalarField& T, const label patchi)...";

	tmp<scalarField> tCv(new scalarField(T.size()));
	scalarField& cv = tCv();

	forAll(T, facei)
	{
		cv[facei] = this->patchFaceMixture(patchi, facei).Cv(T[facei]);
	}

	Info << " end.." << endl;
	return tCv;
}


template<class MixtureType>
Foam::tmp<Foam::volScalarField> Foam::hPdfThermo<MixtureType>::Cv() const
{
	Info << "pdfThermo Cv()...";

	const fvMesh& mesh = T_.mesh();

	tmp<volScalarField> tCv
	(
		new volScalarField
		(
		    IOobject
		    (
			"Cv",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		    ),
		    mesh,
		    dimensionSet(0, 2, -2, -1, 0),
		    T_.boundaryField().types()
		)
	);

	volScalarField& cv = tCv();

	forAll(T_, celli)
	{
		cv[celli] = this->cellMixture(celli).Cv(T_[celli]);
	}

	forAll(T_.boundaryField(), patchi)
	{
		const fvPatchScalarField& pT = T_.boundaryField()[patchi];
		fvPatchScalarField& pCv = cv.boundaryField()[patchi];

		forAll(pT, facei)
		{
		    pCv[facei] = this->patchFaceMixture(patchi, facei).Cv(pT[facei]);
		}
	}

	Info << " end.." << endl;
	return tCv;
}


template<class MixtureType>
bool Foam::hPdfThermo<MixtureType>::read()
{
    if (basicPdfThermo::read())
    {
        MixtureType::read(*this);
        return true;
    }
    else
    {
        return false;
    }
}

template<class MixtureType>
void Foam::hPdfThermo<MixtureType>::update()
{
	std::vector<double> extracted(7);
	double csiv2_normalized = 0.;
	double defect = 0.; 


	const scalarField& csi 		= csi_.internalField();
	const scalarField& csiv2 	= csiv2_.internalField();
	const scalarField& chi_st 	= chi_st_.internalField();
	const scalarField& HCells 	= H_.internalField();

	scalarField& TCells 		= T_.internalField();
	scalarField& RhoCells 		= density_reynolds_.internalField();
	scalarField& asCells 		= as_.internalField();
	scalarField& muCells 		= mu_favre_.internalField();
	scalarField& alphaCells 	= alpha_favre_.internalField();

	double small_eps = 1.e-6;
	double small_chi_st = 1.e-8;
 	
	// Internal cells 
	forAll(csi, celli)
	{
		double max_chi = max(small_chi_st,chi_st[celli]);

		if (adiabaticMode == false)	
			defect = HCells[celli] - (HOxidizer+csi[celli]*(HFuel-HOxidizer));

		if (csi[celli]<=small_eps)		// Pure oxidizer
		{
			flamelets_library.GetMeanValues(0., 0., max_chi, defect, extracted);
		}			

		else if (csi[celli]>=(1.-small_eps))	// Pure fuel	
		{
			flamelets_library.GetMeanValues(1., 0., max_chi, defect, extracted);
		}			

		else					// Mixture
		{
			csiv2_normalized = csiv2[celli] / (csi[celli]*(1.-csi[celli]));	// Normalized mixture fraction variance
	
			if (csiv2_normalized >= 0.98)
				flamelets_library.GetMeanValues(csi[celli], 0.98, max_chi, defect, extracted);
			else if (csiv2_normalized < 0.)
				flamelets_library.GetMeanValues(csi[celli], 0.00, max_chi, defect, extracted);
			else
				flamelets_library.GetMeanValues(csi[celli], csiv2_normalized, max_chi, defect, extracted);	
		}			

		TCells[celli] 		= extracted[1];
		RhoCells[celli]		= extracted[2]; 
		asCells[celli]		= extracted[3]; 
		muCells[celli]		= extracted[4]; 
		alphaCells[celli]	= extracted[5]; 
	}

	// Boundary conditions
	if (adiabaticMode == true)
	{
		forAll(csi_.boundaryField(), patchi)
		{
			const fvPatchScalarField& pcsi 		= csi_.boundaryField()[patchi];
			const fvPatchScalarField& pcsiv2 	= csiv2_.boundaryField()[patchi];
			const fvPatchScalarField& pchi_st 	= chi_st_.boundaryField()[patchi];

			fvPatchScalarField& pt		= T_.boundaryField()[patchi];
			fvPatchScalarField& prho 	= density_reynolds_.boundaryField()[patchi];
			fvPatchScalarField& pas 	= as_.boundaryField()[patchi];
			fvPatchScalarField& pmu 	= mu_favre_.boundaryField()[patchi];
			fvPatchScalarField& palpha 	= alpha_favre_.boundaryField()[patchi];

			forAll(pcsi, facei)
			{

				double max_chi = max(small_chi_st, pchi_st[facei]);

				if (pcsi[facei]<=0.)		// Pure oxidizer
				{
					flamelets_library.GetMeanValues(0., 0., max_chi, defect, extracted);
				}			

				else if (pcsi[facei]>=1.)	// Pure fuel	
				{
					flamelets_library.GetMeanValues(1., 0., max_chi, defect, extracted);
				}			

				else				// Mixture
				{
					csiv2_normalized = pcsiv2[facei] / (pcsi[facei]*(1.-pcsi[facei]));	// Normalized mixture fraction variance

					if (csiv2_normalized >= 0.98)
						flamelets_library.GetMeanValues(pcsi[facei], 0.98, max_chi, defect, extracted);
					else if (csiv2_normalized < 0.)
						flamelets_library.GetMeanValues(pcsi[facei], 0.00, max_chi, defect, extracted);
					else
						flamelets_library.GetMeanValues(pcsi[facei], csiv2_normalized, max_chi,defect, extracted);
				}			

				pt[facei] 	= extracted[1];
				prho[facei]	= extracted[2];
				pas[facei]	= extracted[3];
				pmu[facei]	= extracted[4];
				palpha[facei]	= extracted[5];

			}
		}
	}
	else
	{
		forAll(csi_.boundaryField(), patchi)
		{

			if (patch_type_T[patchi] == 0)
			{
				const fvPatchScalarField& pcsi 		= csi_.boundaryField()[patchi];
				const fvPatchScalarField& pcsiv2 	= csiv2_.boundaryField()[patchi];
				const fvPatchScalarField& pchi_st 	= chi_st_.boundaryField()[patchi];
				const fvPatchScalarField& ph		= H_.boundaryField()[patchi];

				fvPatchScalarField& pt		= T_.boundaryField()[patchi];
				fvPatchScalarField& prho 	= density_reynolds_.boundaryField()[patchi];
				fvPatchScalarField& pas 	= as_.boundaryField()[patchi];
				fvPatchScalarField& pmu 	= mu_favre_.boundaryField()[patchi];
				fvPatchScalarField& palpha 	= alpha_favre_.boundaryField()[patchi];

				forAll(pcsi, facei)
				{

					double max_chi = max(small_chi_st, pchi_st[facei]);

					defect = ph[facei] - (HOxidizer+pcsi[facei]*(HFuel-HOxidizer));

					if (pcsi[facei]<=0.)		// Pure oxidizer
					{
						flamelets_library.GetMeanValues(0., 0., max_chi, defect, extracted);
					}			

					else if (pcsi[facei]>=1.)	// Pure fuel	
					{
						flamelets_library.GetMeanValues(1., 0., max_chi, defect, extracted);
					}			

					else						// Mixture
					{
						csiv2_normalized = pcsiv2[facei] / (pcsi[facei]*(1.-pcsi[facei]));	// Normalized mixture fraction variance

						if (csiv2_normalized >= 0.98)
							flamelets_library.GetMeanValues(pcsi[facei], 0.98, max_chi, defect, extracted);
						else if (csiv2_normalized < 0.)
							flamelets_library.GetMeanValues(pcsi[facei], 0.00, max_chi, defect, extracted);
						else
							flamelets_library.GetMeanValues(pcsi[facei], csiv2_normalized, max_chi,defect, extracted);
					}			

					pt[facei] 	= extracted[1];
					prho[facei]	= extracted[2];
					pas[facei]	= extracted[3];
					pmu[facei]	= extracted[4];
					palpha[facei]	= extracted[5];
				}
			}
			// Added for BC - fixedTemperature
			else if (patch_type_T[patchi] == 1 && patch_type_H[patchi] == 1 && patch_type_csi[patchi] == 1)	// fixed temperature and fixed enhalpy --> INLETS
			{
				const fvPatchScalarField& pcsi 		= csi_.boundaryField()[patchi];
				const fvPatchScalarField& pcsiv2 	= csiv2_.boundaryField()[patchi];
				const fvPatchScalarField& pchi_st 	= chi_st_.boundaryField()[patchi];

				fvPatchScalarField& prho 		= density_reynolds_.boundaryField()[patchi];
				fvPatchScalarField& pas 		= as_.boundaryField()[patchi];
				fvPatchScalarField& pmu 		= mu_favre_.boundaryField()[patchi];
				fvPatchScalarField& palpha 		= alpha_favre_.boundaryField()[patchi];

				forAll(pcsi, facei)
				{

					defect = 0.; // adiabat inlet

					double max_chi = max(small_chi_st, pchi_st[facei]);

					if (pcsi[facei]<=0.)		// Pure oxidizer
					{
						flamelets_library.GetMeanValues(0., 0., max_chi, defect, extracted);
					}			

					else if (pcsi[facei]>=1.)	// Pure fuel	
					{
						flamelets_library.GetMeanValues(1., 0., max_chi, defect, extracted);
					}			

					else				// Mixture
					{
						csiv2_normalized = pcsiv2[facei] / (pcsi[facei]*(1.-pcsi[facei]));	// Normalized mixture fraction variance

						if (csiv2_normalized >= 0.98)
						{
							flamelets_library.GetMeanValues(pcsi[facei], 0.98, max_chi, defect, extracted);
						}					
						else if (csiv2_normalized < 0.)
						{
							flamelets_library.GetMeanValues(pcsi[facei], 0.00, max_chi, defect, extracted);
						}					
						else
						{
							flamelets_library.GetMeanValues(pcsi[facei], csiv2_normalized, max_chi,defect, extracted);
						}
					}			

					prho[facei]	= extracted[2];
					pas[facei]	= extracted[3];
					pmu[facei]	= extracted[4];
					palpha[facei]	= extracted[5];
				}
			}
			// Added for BC - fixedTemperature
			else if (patch_type_T[patchi] == 1 && patch_type_H[patchi] == 1 && patch_type_csi[patchi] == 0)	// fixed temperature and fixedEnthalpie - no inlet
			{

				const fvPatchScalarField& pcsi 		= csi_.boundaryField()[patchi];
				const fvPatchScalarField& pcsiv2 	= csiv2_.boundaryField()[patchi];
				const fvPatchScalarField& pchi_st 	= chi_st_.boundaryField()[patchi];

				fvPatchScalarField& ph			= H_.boundaryField()[patchi];
				fvPatchScalarField& pt			= T_.boundaryField()[patchi];
				fvPatchScalarField& prho 		= density_reynolds_.boundaryField()[patchi];
				fvPatchScalarField& pas 		= as_.boundaryField()[patchi];
				fvPatchScalarField& pmu 		= mu_favre_.boundaryField()[patchi];
				fvPatchScalarField& palpha 		= alpha_favre_.boundaryField()[patchi];

				forAll(pcsi, facei)
				{

					double max_chi = max(small_chi_st, pchi_st[facei]);

					if (pcsi[facei]<=0.)		// Pure oxidizer
					{
						defect = flamelets_library.GetEnthalpyDefectFromTemperature(0., 0., max_chi, pt[facei]);	
						ph[facei] = defect + HOxidizer;
						flamelets_library.GetMeanValues(0., 0., max_chi, defect, extracted);
					}			

					else if (pcsi[facei]>=1.)	// Pure fuel	
					{
						defect = flamelets_library.GetEnthalpyDefectFromTemperature(1., 0., max_chi, pt[facei]);
						ph[facei] = defect + HFuel;
						flamelets_library.GetMeanValues(1., 0., max_chi, defect, extracted);
					}			

					else				// Mixture
					{
						csiv2_normalized = pcsiv2[facei] / (pcsi[facei]*(1.-pcsi[facei]));	// Normalized mixture fraction variance

						if (csiv2_normalized >= 0.98)
						{
							defect = flamelets_library.GetEnthalpyDefectFromTemperature(pcsi[facei], 0.98, max_chi, pt[facei]);	
							ph[facei] = defect + (HOxidizer+pcsi[facei]*(HFuel-HOxidizer));
							flamelets_library.GetMeanValues(pcsi[facei], 0.98, max_chi, defect, extracted);
						}					
						else if (csiv2_normalized < 0.)
						{
							defect = flamelets_library.GetEnthalpyDefectFromTemperature(pcsi[facei], 0.00, max_chi, pt[facei]);	
							ph[facei] = defect + (HOxidizer+pcsi[facei]*(HFuel-HOxidizer));
							flamelets_library.GetMeanValues(pcsi[facei], 0.00, max_chi, defect, extracted);
						}					
						else
						{
							defect = flamelets_library.GetEnthalpyDefectFromTemperature(pcsi[facei], csiv2_normalized, max_chi, pt[facei]);	
							ph[facei] = defect + (HOxidizer+pcsi[facei]*(HFuel-HOxidizer));
							flamelets_library.GetMeanValues(pcsi[facei], csiv2_normalized, max_chi,defect, extracted);
						}
					}			

					prho[facei]	= extracted[2];
					pas[facei]	= extracted[3];
					pmu[facei]	= extracted[4];
					palpha[facei]	= extracted[5];
				} 
			}
			
			// For wrong BC give an error
			else 
			{
	//			Info << "Temperature: " << patch_type_T[patchi] << endl;
	//			Info << "Enthalpie: "   << patch_type_H[patchi] << endl;
	//			Info << "Mischungsbruch: " << patch_type_csi[patchi] << endl; 
			        FatalErrorIn
            			(
			                "hPdfThermo<MixtureType>::update()"
			        )   
				<< "Boundary conditions are wrong: "
				<< "fixed temperature BC must be fixed enthaplie BC; the value will be overwritten from the code"
		                << abort(FatalError);
			} 


		} // patch cycle

	} // non adiabatic mode
}

template<class MixtureType>
void Foam::hPdfThermo<MixtureType>::updateMassFractions()
{
	std::vector<double> extracted(flamelets_library.number_of_species()+1);
	double csiv2_normalized = 0.;
	double defect = 0.; 

	const scalarField& csi 		= csi_.internalField();
	const scalarField& csiv2 	= csiv2_.internalField();
	const scalarField& chi_st 	= chi_st_.internalField();
	const scalarField& HCells 	= H_.internalField();

	double small_eps = 1.e-6;
	double small_chi_st = 1.e-8;

	forAll(csi, celli)
	{
		double max_chi = max(small_chi_st,chi_st[celli]);

		if (adiabaticMode == false)	
		defect = HCells[celli] - (HOxidizer+csi[celli]*(HFuel-HOxidizer));

		if (csi[celli]<=small_eps)			// Pure oxidizer
		{
			flamelets_library.ExtractMeanValues(0., 0., max_chi, defect, extracted);
		}			

		else if (csi[celli]>=(1.-small_eps))	// Pure fuel	
		{
			flamelets_library.ExtractMeanValues(1., 0., max_chi, defect, extracted);
		}			

		else									// Mixture
		{
			csiv2_normalized = csiv2[celli] / (csi[celli]*(1.-csi[celli]));	// Normalized mixture fraction variance

			if (csiv2_normalized >= 0.98)
				flamelets_library.ExtractMeanValues(csi[celli], 0.98, max_chi, defect, extracted);
			else if (csiv2_normalized < 0.)
				flamelets_library.ExtractMeanValues(csi[celli], 0.00, max_chi, defect, extracted);
			else
				flamelets_library.ExtractMeanValues(csi[celli], csiv2_normalized, max_chi, defect, extracted);	
		}			

		for(int j=0;j<flamelets_library.number_of_species()+1;j++)
		{
			if(j<flamelets_library.number_of_species()) 
			{
				omega_[j].internalField()[celli] = extracted[j+1];
			}
/*			else
			{
				omega_[j].internalField()[celli] =     1/(omega_[0].internalField()[celli]/28.0104
									+ omega_[1].internalField()[celli]/2.0158
									+ omega_[2].internalField()[celli]/18.0152
									+ omega_[3].internalField()[celli]/32.9988
									+ omega_[4].internalField()[celli]/28.0134 
									+ omega_[5].internalField()[celli]/44.0098
									+ omega_[6].internalField()[celli]/17.0073)
									*1e6*omega_[0].internalField()[celli]/28.0104;
			}*/
		}
	}

	forAll(csi_.boundaryField(), patchi)
	{
		const fvPatchScalarField& pcsi 		= csi_.boundaryField()[patchi];
		const fvPatchScalarField& pcsiv2 	= csiv2_.boundaryField()[patchi];
		const fvPatchScalarField& pchi_st 	= chi_st_.boundaryField()[patchi];
		const fvPatchScalarField& ph		= H_.boundaryField()[patchi];

		forAll(pcsi, facei)
		{

			double max_chi = max(small_chi_st, pchi_st[facei]);

			if (adiabaticMode == false)	
				defect = ph[facei] - (HOxidizer+pcsi[facei]*(HFuel-HOxidizer));

			if (pcsi[facei]<=small_eps)		// Pure oxidizer
			{
				flamelets_library.ExtractMeanValues(0., 0., max_chi, defect, extracted);
			}			

			else if (pcsi[facei]>=(1.-small_eps))	// Pure fuel	
			{
				flamelets_library.ExtractMeanValues(1., 0., max_chi, defect, extracted);
			}			

			else					// Mixture
			{
				csiv2_normalized = pcsiv2[facei] / (pcsi[facei]*(1.-pcsi[facei]));	// Normalized mixture fraction variance

				if (csiv2_normalized >= 0.98)
					flamelets_library.ExtractMeanValues(pcsi[facei], 0.98, max_chi, defect, extracted);
				else if (csiv2_normalized < 0.)
					flamelets_library.ExtractMeanValues(pcsi[facei], 0.00, max_chi, defect, extracted);
				else
					flamelets_library.ExtractMeanValues(pcsi[facei], csiv2_normalized, max_chi,defect, extracted);	
			}			

			for(int j=0;j<flamelets_library.number_of_species()+1;j++)
			{
				if(j<flamelets_library.number_of_species())
				{
					omega_[j].boundaryField()[patchi][facei] = extracted[j+1];
				}
/*				else
				{
					omega_[j].boundaryField()[patchi][facei] = omega_[0].boundaryField()[patchi][facei]; //1/(omega_[0]/28.0104+omega_[1]/2.0158+omega_[2]/18.0152+omega_[3]/32.9988+omega_[4]/28.0134+omega[5]/44.0098+omega[6]/17.0073)*1e6*omega_[0]/28.0104;
				}*/
			}
		}
	}

}

template<class MixtureType>
void Foam::hPdfThermo<MixtureType>::ErrorMessage(const string message)
{
	Info << "Class: pdfThermo" << endl;
	Info << "Error: " << message << endl;
	getchar();
	//exit(-1);
}


// ************************************************************************* //
