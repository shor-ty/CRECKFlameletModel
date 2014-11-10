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

#ifndef OpenSMOKE_PDF_NonAdiabaticFlamelet_Library_H
#define OpenSMOKE_PDF_NonAdiabaticFlamelet_Library_H

#include "OpenSMOKE_PDF_Flamelet.hpp"
#include "OpenSMOKE_PDF_Flamelet_Library.hpp"
#include "OpenSMOKE_ChiDistribution.hpp"

class OpenSMOKE_PDF_NonAdiabaticFlamelet_Library
{

public:

	OpenSMOKE_PDF_NonAdiabaticFlamelet_Library();
	void Read();
	void Summary();
	
	void   GetMeanValues(const double csi, const double csiv2, const double chi_st, const double phi, std::vector<double>& extracted);
	void   ExtractMeanValues(const double csi, const double csiv2, const double chi_st, const double phi, std::vector<double> &omegaFavre);
	double GetEnthalpyDefectFromTemperature(const double csi, const double csiv2, const double chi_st, const double T);
	
	void SetLibraryPath(const std::string path_library);
	void UnsetExcludeColdFlamelets();
	void SetLogNormalChiDistribution(const double sigma, const int number_of_points);
	void SetAdiabaticMode();
	void SetSpeciesToExtract(const std::vector<std::string> names);
	void SetSpeciesToExtract(const std::string list_of_names);
	int number_of_species();
	const std::vector<std::string> species();
	int index_of_species(const std::string name);

private:

	std::string			_path_library;						// library name
	bool			_iExcludeColdFlamelets;
	chiPDF_kinds	_iChiPDF;

	OpenSMOKE_PDF_Flamelet_Library *flamelet_libraries;

private:

	std::string _name;									// library name
	int	   _nphi;	

	int				number_of_species_to_extract_;	//
	std::vector<std::string>	names_of_species_to_extract_;	//
	std::vector<double>  _phi;

private:

	bool   _iMultiScalarDissipationRates;
	bool   _adiabatic_mode;
	
	double _temperature_f_fuel;
	double _temperature_f_oxidizer;
	double _enthalpy_f_fuel;
	double _enthalpy_f_oxidizer;
	double _density_r_fuel;	
	double _density_r_oxidizer;	
	double _phi_max;
	double _phi_min;
	int    _jPhiAdiabatic;

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

#endif // OpenSMOKE_PDF_Flamelet_NonAdibaticLibrary_H
	
