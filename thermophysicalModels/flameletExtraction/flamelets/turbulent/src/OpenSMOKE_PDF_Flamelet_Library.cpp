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
#include <string.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

#include "OpenSMOKE_PDF_Flamelet_Library.hpp"
#include "OpenSMOKE_PDF_Flamelet.hpp"
#include "OpenSMOKE_ChiDistribution.hpp"

OpenSMOKE_PDF_Flamelet_Library::OpenSMOKE_PDF_Flamelet_Library()
{
	_name = "[not assigned]";
	_path_library = "constant/PDF-Library";
	_temperature_f_max = -1.e16;
	_temperature_f_min =  1.e16;
	_density_r_max = -1.e16;
	_density_r_min =  1.e16;
	_number_of_species_to_extract = 0;
	_iMultiScalarDissipationRates = true;
	_showFlamelet = false;

	_iExcludeColdFlamelets = true;
	_iChiPDF = CHI_PDF_DIRAC;
	_chi_log_normal_sigma = 1.31;
	_chi_log_normal_number_of_points = 40;
}

void OpenSMOKE_PDF_Flamelet_Library::SetLibraryPath(const string path_library)
{
	_path_library = path_library;
}

void OpenSMOKE_PDF_Flamelet_Library::UnsetExcludeColdFlamelets()
{
	_iExcludeColdFlamelets = false;
}

void OpenSMOKE_PDF_Flamelet_Library::SetShowFlamelet()
{
	_showFlamelet = true;
}

void OpenSMOKE_PDF_Flamelet_Library::SetLogNormalChiDistribution(const double sigma, const int number_of_points)
{
	_iChiPDF = CHI_PDF_LOG_NORMAL;
	_chi_log_normal_sigma = sigma;
	_chi_log_normal_number_of_points = number_of_points;
}

void OpenSMOKE_PDF_Flamelet_Library::SetSpeciesToExtract(const vector<string> names)
{
	_number_of_species_to_extract = names.size()-1;
	_names_of_species_to_extract = names;
}

void OpenSMOKE_PDF_Flamelet_Library::Read()
{
	int i;
	string tag;

	_name = _path_library + "/LookUpTable.out";

	ifstream fInput;
	fInput.open(_name.c_str(), ios::in);

	fInput >> tag;
	if (tag != "Adiabatic")
		ErrorMessage("Expected: Adiabatic, Found: " + tag);

	fInput >> tag;
	if (tag != "nc")
		ErrorMessage("Expected: nc, Found: " + tag);
	fInput >> _nc;

	fInput >> tag;
	if (tag != "chi")
		ErrorMessage("Expected: chi, Found: " + tag);
	fInput >> _n;

	_chi_st.resize(_n+1);
	for(i=1;i<=_n;i++)
		fInput >> _chi_st[i];

	fInput.close();
	
	_flame.resize(_n+1);
	for(i=1;i<=_n;i++)
	{
		stringstream index; index << i;
		string flamelet_name = _path_library + "/SR_" + index.str() + ".bin";
		if (_number_of_species_to_extract > 0)
			_flame[i].SetSpeciesToExtract(_names_of_species_to_extract);
		_flame[i].ReadBinary(flamelet_name, i, _chi_st[i]);
		if(_showFlamelet == true)
		{
			_flame[i].Summary();
		}
	}
	
	for(i=1;i<=_n;i++)
	{
		if (_temperature_f_max < _flame[i].temperature_f_max()) 
			_temperature_f_max = _flame[i].temperature_f_max();

			if (_temperature_f_min > _flame[i].temperature_f_min()) 
				_temperature_f_min = _flame[i].temperature_f_min();

			if (_density_r_max < _flame[i].density_r_max()) 
				_density_r_max = _flame[i].density_r_max();

			if (_density_r_min > _flame[i].density_r_min()) 
				_density_r_min = _flame[i].density_r_min();
		}

	for(i=1;i<=_n-1;i++)
		if (_chi_st[i] >= _chi_st[i+1])
			ErrorMessage("Scalar dissipation rates are not in the right order...");
	
	if (_iChiPDF == CHI_PDF_DIRAC && _iExcludeColdFlamelets == true)
	{
		for(i=1;i<=_n;i++)
			if (_flame[i].cold() == true)
			{
				_n=_n-1;
				break;
			}
		_chi_st.resize(_n+1);
		_flame.resize(_n+1);
	}

	_chi_st_min = _chi_st[1];
	_chi_st_max = _chi_st[_n];

 	if (_n==1)	
		_iMultiScalarDissipationRates = false;
	
	if (_n<3 && _iChiPDF == CHI_PDF_LOG_NORMAL)
		ErrorMessage("Log-normal distribution function for scalar dissipation rate can be used only with at least three flamelets...");
		
	// Prepare Scalar Dissipation rate distribution function
	if (_iChiPDF == CHI_PDF_LOG_NORMAL)
	{
		// Scalar Dissipation Rate Probability distribution function options
		_chi_pdf.SetSigma(_chi_log_normal_sigma);
		_chi_pdf.SetNumberOfPoints(_chi_log_normal_number_of_points);
		//_chi_pdf.SetLowerLimit(1.e-6);
		//_chi_pdf.SetHigherLimit(1.e3);
		//_chi_pdf.SetFixedPointRatio();
		//_chi_pdf.SetAccurateCalculation();
		//_chi_pdf.SetWeightThreshold(1.e-10);

		_chi_pdf.BuildGrid();
		_chi_pdf.AssignScalarDissipationRates(_chi_st);
		_t_favre.resize(_n+1);
		_rho_reynolds.resize(_n+1);
		_enthalpy_favre.resize(_n+1);
		_cp_favre.resize(_n+1);
		_mu_favre.resize(_n+1);
		_as_favre.resize(_n+1);
		_alpha_favre.resize(_n+1);
		
		if (_number_of_species_to_extract != 0)
		{
			_w_favre.resize(_n+1);
			for(int j=1;j<=_n;j++)
				_w_favre[j].resize(_number_of_species_to_extract+1);
		}
	}
	
	_temperature_f_fuel 		= _flame[1].temperature_f_fuel();
	_temperature_f_oxidizer 	= _flame[1].temperature_f_oxidizer();
	_enthalpy_f_fuel 			= _flame[1].enthalpy_f_fuel();
	_enthalpy_f_oxidizer 		= _flame[1].enthalpy_f_oxidizer();
	_density_r_fuel 			= _flame[1].density_r_fuel();
	_density_r_oxidizer 		= _flame[1].density_r_oxidizer();
}

void OpenSMOKE_PDF_Flamelet_Library::GetMeanValues(const double csi, const double csiv2, const double chi_st, vector<double>& extracted)
{
	if (_iMultiScalarDissipationRates == true)
	{
		if(_iChiPDF == CHI_PDF_LOG_NORMAL)	
		{	
			GetMeanValuesMultiScalarDissipationRatesLogNormal(csi, csiv2, chi_st, extracted);
//			cout << "Hier" << endl;
		}
		else
			GetMeanValuesMultiScalarDissipationRatesDirac(csi, csiv2, chi_st, extracted);
	}
	else											
		_flame[1].GetMeanValues(csi, csiv2, extracted);
}

void OpenSMOKE_PDF_Flamelet_Library::GetMeanValuesMultiScalarDissipationRatesDirac(const double csi, const double csiv2, const double chi_st, vector<double>& extracted)
{
	if (chi_st <= _chi_st_min)	
	{
		_flame[1].GetMeanValues(csi, csiv2, extracted);
		return;
	}
	else if (chi_st >= _chi_st_max)	
	{
		_flame[_n].GetMeanValues(csi, csiv2, extracted);
		return;
	}
	else
	{
		int kChi = 0;
		double ratio;
		vector<double> extracted_1(extracted.size());
		vector<double> extracted_2(extracted.size());
		
		// Finding kChi
		for(int k=2;k<=_n;k++)
			if (chi_st <= _chi_st[k])
			{
				kChi = k-1;
				break;
			}
			
		// Finding flamelets mean values
		_flame[kChi].GetMeanValues(csi, csiv2, extracted_1);
		_flame[kChi+1].GetMeanValues(csi, csiv2, extracted_2);
		
		// Interpolation ratio
		ratio = (chi_st - _chi_st[kChi])/(_chi_st[kChi+1] - _chi_st[kChi]);
		
		// Mean values
		for(unsigned int k=1;k<=extracted.size()-1;k++)
			extracted[k] = extracted_1[k] + ratio*(extracted_2[k]-extracted_1[k]);
	}
}

void OpenSMOKE_PDF_Flamelet_Library::GetMeanValuesMultiScalarDissipationRatesLogNormal(const double csi, const double csiv2, const double chi_st, vector<double>& extracted)
{
	_chi_pdf.AssignMeanScalarDissipationRate(chi_st);

	for(int j=_chi_pdf.iXinf();j<=_chi_pdf.iXsup();j++)
	{
		_flame[j].GetMeanValues(csi, csiv2, extracted);
		_t_favre[j] 		= extracted[1];
		_rho_reynolds[j] 	= extracted[2];
		_as_favre[j] 		= extracted[3];
		_mu_favre[j] 		= extracted[4];
		_alpha_favre[j] 	= extracted[5];
		_cp_favre[j] 		= extracted[6];
	}

	extracted[1]	= _chi_pdf.GetMeanValue(_t_favre);
	extracted[2]	= _chi_pdf.GetMeanValue(_rho_reynolds);
	extracted[3]	= _chi_pdf.GetMeanValue(_as_favre);
	extracted[4]	= _chi_pdf.GetMeanValue(_mu_favre);
	extracted[5]	= _chi_pdf.GetMeanValue(_alpha_favre);
	extracted[6]	= _chi_pdf.GetMeanValue(_cp_favre);
}


void OpenSMOKE_PDF_Flamelet_Library::ExtractMeanValues(const double csi, const double csiv2, const double chi_st, vector<double> &omegaFavre)
{
	if (_iMultiScalarDissipationRates == true)
	{
		if(_iChiPDF == CHI_PDF_LOG_NORMAL)		
			ExtractMeanValuesMultiScalarDissipationRatesLogNormal(csi, csiv2, chi_st, omegaFavre);
		else
			ExtractMeanValuesMultiScalarDissipationRatesDirac(csi, csiv2, chi_st, omegaFavre);
	}
	else											
		_flame[1].ExtractMeanValues(csi, csiv2, omegaFavre);
}

void OpenSMOKE_PDF_Flamelet_Library::ExtractMeanValuesMultiScalarDissipationRatesDirac(const double csi, const double csiv2, const double chi_st, vector<double> &omegaFavre)
{
	if (chi_st <= _chi_st_min)	
	{
		_flame[1].ExtractMeanValues(csi, csiv2, omegaFavre);
		return;
	}
	else if (chi_st >= _chi_st_max)	
	{
		_flame[_n].ExtractMeanValues(csi, csiv2, omegaFavre);
		return;
	}
	else
	{
		int kChi = 0;
		double ratio;
		vector<double> omegaFavre_1, omegaFavre_2;
		omegaFavre_1.resize(_number_of_species_to_extract+1);
		omegaFavre_2.resize(_number_of_species_to_extract+1);
		
		// Finding kChi
		for(int k=2;k<=_n;k++)
			if (chi_st <= _chi_st[k])
			{
				kChi = k-1;
				break;
			}
			
		// Finding flamelets mean values
		_flame[kChi].ExtractMeanValues(csi, csiv2, omegaFavre_1);
		_flame[kChi+1].ExtractMeanValues(csi, csiv2, omegaFavre_2);
		
		// Interpolation ratio
		ratio = (chi_st - _chi_st[kChi])/(_chi_st[kChi+1] - _chi_st[kChi]);
		
		// Mean values
		for(int j=1; j<=_number_of_species_to_extract;j++)
			omegaFavre[j] = omegaFavre_1[j] + ratio*(omegaFavre_2[j]-omegaFavre_1[j]);
	}
}

void OpenSMOKE_PDF_Flamelet_Library::ExtractMeanValuesMultiScalarDissipationRatesLogNormal(const double csi, const double csiv2, const double chi_st, vector<double> &omegaFavre)
{
	_chi_pdf.AssignMeanScalarDissipationRate(chi_st);

	for(int j=_chi_pdf.iXinf();j<=_chi_pdf.iXsup();j++)
		_flame[j].ExtractMeanValues(csi, csiv2, _w_favre[j]);

	omegaFavre = _chi_pdf.ExtractMeanValue(_w_favre);
}

void OpenSMOKE_PDF_Flamelet_Library::Summary(int i, int max)
{
	cout << endl;
	cout << "Flamelet library property:" << endl << endl;
	cout << "     + File name:                         " << _name << endl;
	cout << "     + Number of Flamelets:               " << _n << endl;
	cout << "     + Minimum stoic. scalar diss. rate:  " << _chi_st_min  << " Hz" << endl;
	cout << "     + Maximum stoic. scalar diss. rate:  " << _chi_st_max  << " Hz" <<  endl;
	cout << "     + Temperature min:                   " << _temperature_f_min << " K" << endl;
	cout << "     + Temperature max:                   " << _temperature_f_max << " K" << endl;
	cout << "     + Density min:                       " << _density_r_min << " kg/m3" << endl;
	cout << "     + Density max:                       " << _density_r_max << " kg/m3" << endl;
	cout << "     + Temperature fuel:                  " << _temperature_f_fuel << " K" << endl;
	cout << "     + Temperature oxidizer:              " << _temperature_f_oxidizer << " K" << endl;
	cout << "     + Density fuel:                      " << _density_r_fuel << " kg/m3" << endl;
	cout << "     + Density oxidizer:                  " << _density_r_oxidizer << " kg/m3" << endl;
	cout << "     + Enthalpy fuel:                     " << _enthalpy_f_fuel << " J/kg" << endl;
	cout << "     + Enthalpy oxidizer:                 " << _enthalpy_f_oxidizer << " J/kg" << endl << endl;
	cout << "/*------------------------------------------------------------------------------*\\" << endl;
	cout << "|--------->     Finished loading Look-Up-Table number: " << i << " out of " << max << "      <---------|" << endl;
	cout << "\\*------------------------------------------------------------------------------*/" << endl;
}

void OpenSMOKE_PDF_Flamelet_Library::ErrorMessage(const string error_message)
{
	cout << "OpenSMOKE_PDF_Flamelet Error" << endl;
	cout << "File name:     " << _name << endl;
	cout << "Error message: " << error_message << endl;
	exit(-1);
}

void OpenSMOKE_PDF_Flamelet_Library::WarningMessage(const string warning_message)
{
	cout << "OpenSMOKE_PDF_Flamelet Warning" << endl;
	cout << "File name:       " << _name << endl;
	cout << "Warning message: " << warning_message << endl;
	cout << "Press enter to continue..." << endl;
	getchar();
}


	
