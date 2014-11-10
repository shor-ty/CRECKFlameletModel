/***************************************************************************
 *   Copyright (C) 2010 by Alberto Cuoci		     				       *
 *   alberto.cuoci@polimi.it                                               *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

#include "OpenSMOKE_PDF_NonAdiabaticFlamelet_Library.hpp"
#include "OpenSMOKE_PDF_Flamelet.hpp"
#include "OpenSMOKE_ChiDistribution.hpp"

OpenSMOKE_PDF_NonAdiabaticFlamelet_Library::OpenSMOKE_PDF_NonAdiabaticFlamelet_Library()
{
	_name = "[not assigned]";
	_path_library = "constant/PDF-Library";
	 number_of_species_to_extract_ = 0;
	_iMultiScalarDissipationRates = true;
	_iExcludeColdFlamelets = true;
	_iChiPDF = CHI_PDF_DIRAC;
	_adiabatic_mode = false;

	_chi_log_normal_sigma = 1.31;
	_chi_log_normal_number_of_points = 40;
}

void OpenSMOKE_PDF_NonAdiabaticFlamelet_Library::SetLibraryPath(const string path_library)
{
	_path_library = path_library;
}

void OpenSMOKE_PDF_NonAdiabaticFlamelet_Library::UnsetExcludeColdFlamelets()
{
	_iExcludeColdFlamelets = false;
}

void OpenSMOKE_PDF_NonAdiabaticFlamelet_Library::SetLogNormalChiDistribution(const double sigma, const int number_of_points)
{
	_iChiPDF = CHI_PDF_LOG_NORMAL;
	_chi_log_normal_sigma = sigma;
	_chi_log_normal_number_of_points = number_of_points;
}

void OpenSMOKE_PDF_NonAdiabaticFlamelet_Library::SetSpeciesToExtract(const vector<string> names)
{
	number_of_species_to_extract_ = names.size()-1;
	names_of_species_to_extract_ = names;
}

void OpenSMOKE_PDF_NonAdiabaticFlamelet_Library::SetSpeciesToExtract(const string list_of_names)
{
	vector<string> names;
	names.push_back("list");
	stringstream stream(list_of_names);
	while (!stream.eof())
	{
		string dummy;
		stream >> dummy;
		if (dummy == "" || dummy == "\n")
			continue;
		names.push_back(dummy);
	}
	SetSpeciesToExtract(names);
}

void OpenSMOKE_PDF_NonAdiabaticFlamelet_Library::SetAdiabaticMode()
{
	_adiabatic_mode = true;
}

int OpenSMOKE_PDF_NonAdiabaticFlamelet_Library::number_of_species()
{
	return number_of_species_to_extract_;
}

const vector<string> OpenSMOKE_PDF_NonAdiabaticFlamelet_Library::species()
{
	return names_of_species_to_extract_;
}

int OpenSMOKE_PDF_NonAdiabaticFlamelet_Library::index_of_species(const string name)
{
	for(int j=1;j<=number_of_species_to_extract_;j++)
		if (names_of_species_to_extract_[j] == name)
			return j;
	ErrorMessage("Species " + name + " cannot be extracted from the flamelet library");
	return 0;
}

void OpenSMOKE_PDF_NonAdiabaticFlamelet_Library::Read()
{
	int i;
	string tag;

	_name = _path_library + "/LookUpTable.out";

	ifstream fInput;
	fInput.open(_name.c_str(), ios::in);

	fInput >> tag;
	if (tag != "EnthalpyDefects")
		ErrorMessage("Expected: EnthalpyDefects, Found: " + tag);

	fInput >> tag;
	if (tag != "nphi")
		ErrorMessage("Expected: nphi, Found: " + tag);
	fInput >> _nphi;

	vector<string> _phi_names(_nphi+1);
	_phi.resize(_nphi+1);
	for(i=1;i<=_nphi;i++)
	{
		fInput >> _phi[i];	
		fInput >> _phi_names[i];
	}

	fInput.close();

	for(i=1;i<=_nphi-1;i++)
		if (_phi[i] <= _phi[i+1])
			ErrorMessage("Enthalpy defects must be specified in decreasing order!");

	_jPhiAdiabatic = 0;
	for(i=1;i<=_nphi;i++)
		if (fabs(_phi[i]) <= 1.e-6)
		{
			_jPhiAdiabatic = i;
			break;
		}
	if (_jPhiAdiabatic == 0)
		ErrorMessage("Missing adiabatic enthalppy defect!");

	// Maximum enthalpy defect
	_phi_max = _phi[_nphi];

	// Minimum enthalpy defect
	_phi_min = _phi[1];

	// Settings
	if (_nphi == 1)	SetAdiabaticMode();

	flamelet_libraries = new OpenSMOKE_PDF_Flamelet_Library[_nphi+1];

	for(i=1;i<=_nphi;i++)
	{
		string path = _path_library + "/" + _phi_names[i];
		flamelet_libraries[i].SetLibraryPath(path);

		if (_iExcludeColdFlamelets == false)
			flamelet_libraries[i].UnsetExcludeColdFlamelets();

		if (_iChiPDF == CHI_PDF_LOG_NORMAL)
			flamelet_libraries[i].SetLogNormalChiDistribution(_chi_log_normal_sigma, _chi_log_normal_number_of_points);

		if (number_of_species_to_extract_ != 0)
			flamelet_libraries[i].SetSpeciesToExtract(names_of_species_to_extract_);

		flamelet_libraries[i].Read();

		if (_adiabatic_mode == true)
			break;
	}

	_temperature_f_fuel 	= flamelet_libraries[_jPhiAdiabatic].temperature_f_fuel();
	_temperature_f_oxidizer = flamelet_libraries[_jPhiAdiabatic].temperature_f_oxidizer();
	_enthalpy_f_fuel 	= flamelet_libraries[_jPhiAdiabatic].enthalpy_f_fuel();
	_enthalpy_f_oxidizer 	= flamelet_libraries[_jPhiAdiabatic].enthalpy_f_oxidizer();
	_density_r_fuel 	= flamelet_libraries[_jPhiAdiabatic].density_r_fuel();
	_density_r_oxidizer 	= flamelet_libraries[_jPhiAdiabatic].density_r_oxidizer();
}

void OpenSMOKE_PDF_NonAdiabaticFlamelet_Library::GetMeanValues(const double csi, const double csiv2, const double chi_st, const double phi, vector<double>& extracted)
{
	if (_adiabatic_mode == true)
	{
		flamelet_libraries[1].GetMeanValues(csi, csiv2, chi_st, extracted);
	}
	else
	{
		if (phi >= _phi_min)	
		{
			flamelet_libraries[1].GetMeanValues(csi, csiv2, chi_st, extracted);
			return;
		}
		else if (phi <= _phi_max)	
		{
			flamelet_libraries[_nphi].GetMeanValues(csi, csiv2, chi_st, extracted);
			return;
		}
		else
		{
			int k = 0;
			double ratio;
			vector<double> extracted_1(extracted.size());
			vector<double> extracted_2(extracted.size());
			
			// Finding k
			for(int j=2;j<=_nphi;j++)
				if (phi >= _phi[j])
				{
					k = j;
					break;
				}
							
			// Finding flamelets mean values
			flamelet_libraries[k].GetMeanValues(csi, csiv2, chi_st, extracted_1);
			flamelet_libraries[k-1].GetMeanValues(csi, csiv2, chi_st, extracted_2);
			
			// Interpolation ratio
			ratio = (phi - _phi[k])/(_phi[k-1]-_phi[k]);
			
			// Mean values
			for(unsigned int i=1;i<=extracted.size()-1;i++)
				extracted[i] = extracted_1[i] + ratio*(extracted_2[i]-extracted_1[i]);
		}
	}
}

double OpenSMOKE_PDF_NonAdiabaticFlamelet_Library::GetEnthalpyDefectFromTemperature(const double csi, const double csiv2, const double chi_st, const double T)
{
	const int iTemperature = 1;
	vector<double> extracted(7);

//	cout << "List" << endl;
//	for(int j=1;j<=_nphi;j++)
//	{
//		GetMeanValues(csi, csiv2, chi_st, _phi[j], extracted);
//		cout << j << " " << extracted[1] << endl;
//	}
//	cout << "List" << endl;

	// Case: T>T0
	GetMeanValues(csi, csiv2, chi_st, _phi_min, extracted);
	double T0 = extracted[iTemperature];
	if (T>T0) return _phi_min;

	// Case: T<Tmin
	GetMeanValues(csi, csiv2, chi_st, _phi_max, extracted);
	double Tmin = extracted[iTemperature];
	if (T<Tmin) return _phi_max;

	// General case
	{
		double phiA = 0.;
		double phiB = 0.;
		double phiC, TC;

		// Finding k
		for(int j=2;j<=_nphi;j++)
		{
			GetMeanValues(csi, csiv2, chi_st, _phi[j], extracted);
			double Tj = extracted[iTemperature];
			if (Tj <= T)	
			{ 
				phiA = _phi[j]; 
				phiB = _phi[j-1];
				
			//	cout << "R " << " " << " " << T << endl;
			//	cout << "0 " << " " << " " << T0 << endl;
			//	cout << "M " << " " << " " << Tmin << endl;
			//	cout << "A " << j   << " " << phiA << endl;
			//	cout << "B " << j-1 << " " << phiB << endl;
				break;
			}
		}
		
		int kmax = 15;	
		for(int k=1;k<=kmax;k++)
		{
			phiC = 0.50*(phiA+phiB);
			GetMeanValues(csi, csiv2, chi_st, phiC, extracted);
			TC = extracted[iTemperature];

		//	cout << "k" << k << " " << phiC << " " << TC << endl;

			if (T<=TC)
			{
				phiB = phiC;
			//	TB 	= TC;
			}
			else
			{
				phiA 	= phiC;
			//	TA	= TC;
			}

			if ( fabs(TC-T)/T < 1e-5)
				break;
		}
	//	getchar();
		return phiC;
	}
}

void OpenSMOKE_PDF_NonAdiabaticFlamelet_Library::ExtractMeanValues(const double csi, const double csiv2, const double chi_st, const double phi, vector<double> &omegaFavre)
{
	if (_adiabatic_mode == true)
	{
		flamelet_libraries[1].ExtractMeanValues(csi, csiv2, chi_st, omegaFavre);
	}
	else
	{
		if (phi >= _phi_min)	
		{
			flamelet_libraries[1].ExtractMeanValues(csi, csiv2, chi_st, omegaFavre);
			return;
		}
		else if (phi <= _phi_max)	
		{
			flamelet_libraries[_nphi].ExtractMeanValues(csi, csiv2, chi_st, omegaFavre);
			return;
		}
		else
		{
			int k = 0;
			double ratio;
			vector<double> omegaFavre_1(omegaFavre.size());
			vector<double> omegaFavre_2(omegaFavre.size());
			
			// Finding k
			for(int j=2;j<=_nphi;j++)
				if (phi >= _phi[j])
				{
					k = j;
					break;
				}
							
			// Finding flamelets mean values
			flamelet_libraries[k].ExtractMeanValues(csi, csiv2, chi_st, omegaFavre_1);
			flamelet_libraries[k-1].ExtractMeanValues(csi, csiv2, chi_st, omegaFavre_2);
			
			// Interpolation ratio
			ratio = (phi - _phi[k])/(_phi[k-1]-_phi[k]);
			
			// Mean values
			for(int w=1;w<=number_of_species_to_extract_;w++)
				omegaFavre[w] = omegaFavre_1[w] + ratio*(omegaFavre_2[w]-omegaFavre_1[w]);
		}
	}
}


void OpenSMOKE_PDF_NonAdiabaticFlamelet_Library::Summary()
{
	cout << endl;
	cout << "OpenSMOKE_PDF_NonAdiabaticFlamelet_Library" << endl;
	cout << "File name:                         " << _name << endl;
	cout << "Number of Enthalpy defects:        " << _nphi << endl;
	cout << "Minimum enthalpy defect:           " << _phi_min/1000.  << " kJ/kg" << endl;
	cout << "Maximum enthalpy defect:           " << _phi_max/1000.  << " kJ/kg" << endl;
	cout << "Adiabatic flamelets:               " << _jPhiAdiabatic              << endl;
	cout << "Temperature fuel:                  " << _temperature_f_fuel << " K" << endl;
	cout << "Temperature oxidizer:              " << _temperature_f_oxidizer << " K" << endl;
	cout << "Density fuel:                      " << _density_r_fuel << " kg/m3" << endl;
	cout << "Density oxidizer:                  " << _density_r_oxidizer << " kg/m3" << endl;	
	cout << "Enthalpy fuel:                     " << _enthalpy_f_fuel << " J/kg" << endl;
	cout << "Enthalpy oxidizer:                 " << _enthalpy_f_oxidizer << " J/kg" << endl;	
	cout << endl;
}

void OpenSMOKE_PDF_NonAdiabaticFlamelet_Library::ErrorMessage(const string error_message)
{
	cout << "OpenSMOKE_PDF_NonAdiabaticFlamelet_Library Error" << endl;
	cout << "File name:     " << _name << endl;
	cout << "Error message: " << error_message << endl;
	exit(-1);
}

void OpenSMOKE_PDF_NonAdiabaticFlamelet_Library::WarningMessage(const string warning_message)
{
	cout << "OpenSMOKE_PDF_NonAdiabaticFlamelet_Library Warning" << endl;
	cout << "File name:       " << _name << endl;
	cout << "Warning message: " << warning_message << endl;
	cout << "Press enter to continue..." << endl;
	getchar();
}


	
