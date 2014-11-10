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

#ifndef OpenSMOKE_PDF_Flamelet_Library_H
#define OpenSMOKE_PDF_Flamelet_Library_H

#include "OpenSMOKE_PDF_Flamelet.hpp"
#include "OpenSMOKE_ChiDistribution.hpp"
#include "OpenSMOKE_External_Functions.hpp"

class OpenSMOKE_PDF_Flamelet_Library
{

public:

	OpenSMOKE_PDF_Flamelet_Library();
	void Read();
	void Summary();
	
	void GetMeanValues(const double csi, const double csiv2, const double chi_st, std::vector<double>& extracted);
	void ExtractMeanValues(const double csi, const double csiv2, const double chi_st, std::vector<double> &omegaFavre);

	void SetLibraryPath(const std::string path_library);
	void UnsetExcludeColdFlamelets();
	void SetLogNormalChiDistribution(const double sigma, const int number_of_points);
	void SetSpeciesToExtract(const std::vector<std::string> names);

private:

	std::string			 _path_library;						// library name
	bool			 _iExcludeColdFlamelets;
	chiPDF_kinds	 _iChiPDF;

private:

	std::string _name;									// library name
	int    _n;										// number of flamelets
	int    _nc;										// number of species

	int				_number_of_species_to_extract;	//
	std::vector<std::string>	_names_of_species_to_extract;	//


	std::vector<double> _chi_st;						// stoichiometric scalar dissipation rate [Hz]
	std::vector< OpenSMOKE_PDF_Flamelet > _flame;	// Flamelets
	OpenSMOKE_ChiDistribution	 _chi_pdf;	// Scalar dissipation rate probability distribution function

	void GetMeanValuesMultiScalarDissipationRatesDirac(const double csi, const double csiv2, const double chi_st, std::vector<double> &extracted);
	void GetMeanValuesMultiScalarDissipationRatesLogNormal(const double csi, const double csiv2, const double chi_st, std::vector<double> &extracted);

	void ExtractMeanValuesMultiScalarDissipationRatesDirac(const double csi, const double csiv2, const double chi_st, std::vector<double> &omegaFavre);
	void ExtractMeanValuesMultiScalarDissipationRatesLogNormal(const double csi, const double csiv2, const double chi_st, std::vector<double> &omegaFavre);

private:

	bool   _iMultiScalarDissipationRates;

	double _temperature_f_max;
	double _temperature_f_min;
	double _density_r_max;
	double _density_r_min;
	
	double _chi_st_min;
	double _chi_st_max;

	std::vector<double>				_t_favre;
	std::vector<double>				_rho_reynolds;
	std::vector< std::vector<double> >	_w_favre;
	std::vector<double>				_enthalpy_favre;
	std::vector<double>				_cp_favre;
	std::vector<double>				_mu_favre;
	std::vector<double>				_as_favre;
	std::vector<double>				_alpha_favre;
	
		
	double _temperature_f_fuel;
	double _temperature_f_oxidizer;
	double _enthalpy_f_fuel;
	double _enthalpy_f_oxidizer;
	double _density_r_fuel;	
	double _density_r_oxidizer;	

	double 	_chi_log_normal_sigma;
	int 	_chi_log_normal_number_of_points;

private:

	void ErrorMessage(const std::string error_message);
	void WarningMessage(const std::string warning_message);
	
public:

	inline double density_r_fuel()
	{
		return _density_r_fuel;
	};

	inline double density_r_oxidizer()
	{
		return _density_r_oxidizer;
	};
	
	inline double enthalpy_f_fuel()
	{
		return _enthalpy_f_fuel;
	};

	inline double enthalpy_f_oxidizer()
	{
		return _enthalpy_f_oxidizer;
	};
	
	inline double temperature_f_fuel()
	{
		return _temperature_f_fuel;
	};
	
	inline double temperature_f_oxidizer()
	{
		return _temperature_f_oxidizer;
	};	
};

#endif // OpenSMOKE_PDF_Flamelet_Library_H
	
