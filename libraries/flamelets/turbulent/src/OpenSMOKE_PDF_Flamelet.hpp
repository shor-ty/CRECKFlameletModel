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
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#ifndef OpenSMOKE_PDF_Flamelet_H
#define OpenSMOKE_PDF_Flamelet_H

enum extraction_modes {OPENSMOKE_FLAMELET_EXTRACTION_ADIABATIC, OPENSMOKE_FLAMELET_EXTRACTION_NONADIABATIC};

class OpenSMOKE_PDF_Flamelet
{

public:
	
	OpenSMOKE_PDF_Flamelet();
	void Read(const std::string file_name, const int index, const double chi_st);
	void ReadBinary(const std::string file_name, const int index, const double chi_st);
	void Summary();
	
	void GetMeanValues(const double csi, const double csiv2, std::vector<double>& extracted);
	void ExtractMeanValues(const double csi, const double csiv2, std::vector<double> &omegaFavre);

private:

	std::string _name;
	int _index;

	double _chi_st;				// stoichiometric scalar dissipation rate [Hz]
	double _chi_max;			// maximum scalar dissipation rate [Hz]
	double _as;					// nominal strain rate [Hz]

	int _n_csi;					// number of mixture fraction points
	int _n_variance;			// number of mixture fraction variance points;

	std::vector<double>	_csi;						// mixture fraction points (favre)
	std::vector<double>	_variance_normal;			// normal variance points (favre)

	std::vector< std::vector<double> > _mf_r;				// mixture fraction points (reynolds)
	std::vector< std::vector<double> > _mfv_r;			// normal variance points (reynolds)

	std::vector< std::vector<double> > _density_r;		// Reynolds density [kg/m3]
	std::vector< std::vector<double> > _temperature_f;	// Favre temperature [K]

	std::vector< std::vector<double> > _enthalpy_f;	// Favre enthalpy [J/kg]
	std::vector< std::vector<double> > _mw_f;			// Favre temperature [kg/kmol]
	std::vector< std::vector<double> > _cp_f;			// Favre temperature [J/kg]
	std::vector< std::vector<double> > _lambda_f;		// Favre thermal conductivity [J/s/m/K]
	std::vector< std::vector<double> > _mu_f;			// Favre dynamic viscosity [kg/m/s]
	std::vector< std::vector<double> > _as_f;			// Favre absorption coefficient [1/m]
	std::vector< std::vector<double> > _alpha_f;		// Favre thermal diffusivity [m2/s]

	std::vector< std::vector< std::vector<double> > > _w_f;	// Favre mass fractions

	int		_number_of_species;
	int		_number_of_species_to_extract;
	std::vector<std::string>	_names_of_species_to_extract;
	std::vector<double>	_mw_species_to_extract;


private:

	double _temperature_f_max;
	double _temperature_f_min;
	double _density_r_max;
	double _density_r_min;
	bool   _cold;
			
private:

	void ErrorMessage(const std::string error_message);
	void WarningMessage(const std::string warning_message);

public:

	void SetSpeciesToExtract(const std::vector<std::string> names);

	double temperature_f_max();
	double temperature_f_min();
	double density_r_max();
	double density_r_min();
	bool cold();
	double density_r_fuel();
	double density_r_oxidizer();
	double enthalpy_f_fuel();
	double enthalpy_f_oxidizer();
	double temperature_f_fuel();
	double temperature_f_oxidizer();

};

#endif // #ifndef OpenSMOKE_PDF_Flamelet_H

	
