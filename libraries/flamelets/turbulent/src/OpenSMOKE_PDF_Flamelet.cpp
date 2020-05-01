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

#include "OpenSMOKE_PDF_Flamelet.hpp"

OpenSMOKE_PDF_Flamelet::OpenSMOKE_PDF_Flamelet()
{
	_name = "NULL";
	_index = 0;
	_temperature_f_max = -1.e16;
	_temperature_f_min =  1.e16;
	_density_r_max = -1.e16;
	_density_r_min =  1.e16;
	_number_of_species_to_extract = 0;

	_cold = false;
}

void OpenSMOKE_PDF_Flamelet::SetSpeciesToExtract(const vector<string> names)
{
	_number_of_species_to_extract = names.size()-1;
	_names_of_species_to_extract = names;
}

void OpenSMOKE_PDF_Flamelet::Read(const string file_name, const int index, const double chi_st)
{
	int i, j;
	string tag;

	_name = file_name;
	_index = index;
	_chi_st = chi_st;

	ifstream fInput;
	fInput.open(file_name.c_str(), ios::in);

	// 1  - mixture fraction (favre): mf
	// 2  - variance of mixture fraction (favre): mfv

	// 3  - mixture fraction (reynolds): mf_reynolds
	// 4  - variance of mixture fraction (reynolds): mfv_reynolds

	// 4  - temperature (favre) (J/kg): T
	// 5  - density (reynolds) (kg/m3): rho
	
	// 7  - specific enthalpy (favre) (J/kg): h
	// 8  - molecular weight (favre) (kg/kmol): mw
	// 9  - specific heat (favre) (J/kg/K): cp
	// 10 - thermal conductivity (favre) (J/s/m/K): lambda
	// 11 - dynamic viscosity (favre) (kg/m/s): mu
	// 12 - absorption coefficient (favre) (1/m): as 

	// 13 - species mass fractions

	// Read flamelet version
	fInput >> tag;
	fInput >> tag;
	if (tag != "0.1")
		ErrorMessage("Only version 0.1 is supported...");

	// Read total number of species in the flamelet
	fInput >> tag;
	if (tag != "nc")
		ErrorMessage("Expected: nc, Found: " + tag);
	fInput >> _number_of_species;

	fInput >> tag;	// 1
	if (tag != "mf")
		ErrorMessage("Expected: mf, Found: " + tag);

	fInput >> _n_csi;
	_csi.resize(_n_csi+1);
	for(i=1;i<=_n_csi;i++)
		fInput >> _csi[i];

	fInput >> tag;	// 2
	if (tag != "mfv")
		ErrorMessage("Expected: mfv, Found: " + tag);

	fInput >> _n_variance;
	_variance_normal.resize(_n_variance+1);
	for(j=1;j<=_n_variance;j++)
		fInput >> _variance_normal[j];


	// Memory allocation
	_mf_r.resize(_n_variance+1);
	_mfv_r.resize(_n_variance+1);
	_temperature_f.resize(_n_variance+1);
	_density_r.resize(_n_variance+1);
	_enthalpy_f.resize(_n_variance+1);
	_mw_f.resize(_n_variance+1);
	_cp_f.resize(_n_variance+1);
	_lambda_f.resize(_n_variance+1);
	_mu_f.resize(_n_variance+1);
	_as_f.resize(_n_variance+1);
	_alpha_f.resize(_n_variance+1);

	for(j=1;j<=_n_variance;j++)
	{
		_mf_r[j].resize(_n_csi+1);
		_mfv_r[j].resize(_n_csi+1);
		_temperature_f[j].resize(_n_csi+1);
		_density_r[j].resize(_n_csi+1);
		_enthalpy_f[j].resize(_n_csi+1);
		_mw_f[j].resize(_n_csi+1);
		_cp_f[j].resize(_n_csi+1);
		_lambda_f[j].resize(_n_csi+1);
		_mu_f[j].resize(_n_csi+1);
		_as_f[j].resize(_n_csi+1);
		_alpha_f[j].resize(_n_csi+1);		
	}

	fInput >> tag;	// 3
	if (tag != "mf_reynolds")
		ErrorMessage("Expected: mf_reynolds, Found: " + tag);
	for(j=1;j<=_n_variance;j++)
		for(i=1;i<=_n_csi;i++)
			fInput >> _mf_r[j][i];

	fInput >> tag; // 4
	if (tag != "mfv_reynolds")
		ErrorMessage("Expected: mfv_reynolds, Found: " + tag);
	for(j=1;j<=_n_variance;j++)
		for(i=1;i<=_n_csi;i++)
			fInput >> _mfv_r[j][i];

	fInput >> tag;	// 5
	if (tag != "temperature")
		ErrorMessage("Expected: temperature, Found: " + tag);
	for(j=1;j<=_n_variance;j++)
		for(i=1;i<=_n_csi;i++)
			fInput >> _temperature_f[j][i];

	fInput >> tag; // 6
	if (tag != "density")
		ErrorMessage("Expected: density, Found: " + tag);
	for(j=1;j<=_n_variance;j++)
		for(i=1;i<=_n_csi;i++)
			fInput >> _density_r[j][i];

	fInput >> tag;	// 7
	if (tag != "enthalpy")
		ErrorMessage("Expected: enthalpy, Found: " + tag);
	for(j=1;j<=_n_variance;j++)
		for(i=1;i<=_n_csi;i++)
			fInput >> _enthalpy_f[j][i];

	fInput >> tag;	// 8
	if (tag != "mw")
		ErrorMessage("Expected: mw, Found: " + tag);
	for(j=1;j<=_n_variance;j++)
		for(i=1;i<=_n_csi;i++)
			fInput >> _mw_f[j][i];

	fInput >> tag;	// 9
	if (tag != "cp")
		ErrorMessage("Expected: cp, Found: " + tag);
	for(j=1;j<=_n_variance;j++)
		for(i=1;i<=_n_csi;i++)
			fInput >> _cp_f[j][i];

	fInput >> tag;	// 10
	if (tag != "lambda")
		ErrorMessage("Expected: lambda, Found: " + tag);
	for(j=1;j<=_n_variance;j++)
		for(i=1;i<=_n_csi;i++)
			fInput >> _lambda_f[j][i];

	fInput >> tag;	// 11
	if (tag != "mu")
		ErrorMessage("Expected: mu, Found: " + tag);
	for(j=1;j<=_n_variance;j++)
		for(i=1;i<=_n_csi;i++)
			fInput >> _mu_f[j][i];

	fInput >> tag;	// 12
	if (tag != "as")
		ErrorMessage("Expected: as, Found: " + tag);
	for(j=1;j<=_n_variance;j++)
		for(i=1;i<=_n_csi;i++)
			fInput >> _as_f[j][i];
			
	// alpha
	for(j=1;j<=_n_variance;j++)
		for(i=1;i<=_n_csi;i++)
			_alpha_f[j][i] = _lambda_f[j][i]/_density_r[j][i]/_cp_f[j][i];		

	if(_number_of_species_to_extract > 0)
	{
		_mw_species_to_extract.resize(_number_of_species_to_extract+1);
		_w_f.resize(_number_of_species_to_extract+1);
		for(int k=1;k<=_number_of_species_to_extract;k++)
		{
			_mw_species_to_extract[k] = 0.;
			_w_f[k].resize(_n_variance+1);
			for(j=1;j<=_n_variance;j++)
				_w_f[k][j].resize(_n_csi+1);
		}
		
		int _count_species_found = 0;
		for(int k=1;k<=_number_of_species;k++)
		{
			double mw;
			double dummy;

			fInput >> tag;
			fInput >> mw;

			bool wFound = false;
			for(int w=1;w<=_number_of_species_to_extract;w++)
				if (_names_of_species_to_extract[w] == tag)
				{
					for(j=1;j<=_n_variance;j++)
						for(i=1;i<=_n_csi;i++)
							fInput >> _w_f[w][j][i];
					_mw_species_to_extract[w] = mw;

					_count_species_found++;
					wFound = true;
					break;
				}

			if (_count_species_found == _number_of_species_to_extract)
				break;
			
			if (wFound == false)				// TODO
				for(j=1;j<=_n_variance;j++)
						for(i=1;i<=_n_csi;i++)
							fInput >> dummy;
		}
		
		for(int w=1;w<=_number_of_species_to_extract;w++)
			if (_mw_species_to_extract[w] == 0.)
				ErrorMessage("The following species is not available in the flamelet library: " + _names_of_species_to_extract[w]);
	}

	fInput.close();


	for(j=1;j<=_n_variance;j++)
		for(i=1;i<=_n_csi;i++)
		{
			if (_temperature_f_max < _temperature_f[j][i]) 
				_temperature_f_max = _temperature_f[j][i];

			if (_temperature_f_min > _temperature_f[j][i]) 
				_temperature_f_min = _temperature_f[j][i];

			if (_density_r_max < _density_r[j][i]) 
				_density_r_max = _density_r[j][i];

			if (_density_r_min > _density_r[j][i]) 
				_density_r_min = _density_r[j][i];
		}

	if ( _temperature_f_max <= 1.05*max(_temperature_f[1][1],_temperature_f[1][_n_csi]) )
		_cold = true;
}

void OpenSMOKE_PDF_Flamelet::ReadBinary(const string file_name, const int index, const double chi_st)
{
	int i, j;
	const int SIZE = 40;
	char tag[SIZE];

	_name = file_name;
	_index = index;
	_chi_st = chi_st;

	fstream fInput;
	fInput.open(file_name.c_str(), ios::in | ios::binary);

	// 1  - mixture fraction (favre): mf
	// 2  - variance of mixture fraction (favre): mfv

	// 3  - mixture fraction (reynolds): mf_reynolds
	// 4  - variance of mixture fraction (reynolds): mfv_reynolds

	// 4  - temperature (favre) (J/kg): T
	// 5  - density (reynolds) (kg/m3): rho
	
	// 7  - specific enthalpy (favre) (J/kg): h
	// 8  - molecular weight (favre) (kg/kmol): mw
	// 9  - specific heat (favre) (J/kg/K): cp
	// 10 - thermal conductivity (favre) (J/s/m/K): lambda
	// 11 - dynamic viscosity (favre) (kg/m/s): mu
	// 12 - absorption coefficient (favre) (1/m): as 

	// 13 - species mass fractions

	// Read flamelet version
	fInput.read(tag, SIZE);
	fInput.read(tag, SIZE);
	if (string(tag) != "0.1")
		ErrorMessage("Only version 0.1 is supported...");

	// Read total number of species in the flamelet
	fInput.read(tag, SIZE);
	if (string(tag) != "nc")
		ErrorMessage("Expected: nc, Found: " + string(tag));
	fInput.read(reinterpret_cast < char * > (&_number_of_species),sizeof(int));

	fInput.read(tag, SIZE);	// 1
	if (string(tag) != "mf")
		ErrorMessage("Expected: mf, Found: " + string(tag));

	fInput.read(reinterpret_cast < char * > (&_n_csi),sizeof(int));
	_csi.resize(_n_csi+1);
	for(i=1;i<=_n_csi;i++)
		fInput.read(reinterpret_cast < char * > (&_csi[i]),sizeof(double));;

	fInput.read(tag, SIZE);	// 2
	if (string(tag) != "mfv")
		ErrorMessage("Expected: mfv, Found: " + string(tag));

	fInput.read(reinterpret_cast < char * > (&_n_variance),sizeof(int));
	_variance_normal.resize(_n_variance+1);
	for(j=1;j<=_n_variance;j++)
		fInput.read(reinterpret_cast < char * > (&_variance_normal[j]),sizeof(double));


	// Memory allocation
	_mf_r.resize(_n_variance+1);
	_mfv_r.resize(_n_variance+1);
	_temperature_f.resize(_n_variance+1);
	_density_r.resize(_n_variance+1);
	_enthalpy_f.resize(_n_variance+1);
	_mw_f.resize(_n_variance+1);
	_cp_f.resize(_n_variance+1);
	_lambda_f.resize(_n_variance+1);
	_mu_f.resize(_n_variance+1);
	_as_f.resize(_n_variance+1);
	_alpha_f.resize(_n_variance+1);

	for(j=1;j<=_n_variance;j++)
	{
		_mf_r[j].resize(_n_csi+1);
		_mfv_r[j].resize(_n_csi+1);
		_temperature_f[j].resize(_n_csi+1);
		_density_r[j].resize(_n_csi+1);
		_enthalpy_f[j].resize(_n_csi+1);
		_mw_f[j].resize(_n_csi+1);
		_cp_f[j].resize(_n_csi+1);
		_lambda_f[j].resize(_n_csi+1);
		_mu_f[j].resize(_n_csi+1);
		_as_f[j].resize(_n_csi+1);
		_alpha_f[j].resize(_n_csi+1);		
	}

	fInput.read(tag, SIZE);	// 3
	if (string(tag) != "mf_reynolds")
		ErrorMessage("Expected: mf_reynolds, Found: " + string(tag));
	for(j=1;j<=_n_variance;j++)
		for(i=1;i<=_n_csi;i++)
			fInput.read(reinterpret_cast < char * > (&_mf_r[j][i]),sizeof(double));

	fInput.read(tag, SIZE);	// 4
	if (string(tag) != "mfv_reynolds")
		ErrorMessage("Expected: mfv_reynolds, Found: " + string(tag));
	for(j=1;j<=_n_variance;j++)
		for(i=1;i<=_n_csi;i++)
			fInput.read(reinterpret_cast < char * > (&_mfv_r[j][i]),sizeof(double));

	fInput.read(tag, SIZE);	// 5
	if (string(tag) != "temperature")
		ErrorMessage("Expected: temperature, Found: " + string(tag));
	for(j=1;j<=_n_variance;j++)
		for(i=1;i<=_n_csi;i++)
			fInput.read(reinterpret_cast < char * > (&_temperature_f[j][i]),sizeof(double));

	fInput.read(tag, SIZE);	// 6
	if (string(tag) != "density")
		ErrorMessage("Expected: density, Found: " + string(tag));
	for(j=1;j<=_n_variance;j++)
		for(i=1;i<=_n_csi;i++)
			fInput.read(reinterpret_cast < char * > (&_density_r[j][i]),sizeof(double));

	fInput.read(tag, SIZE);	// 7
	if (string(tag) != "enthalpy")
		ErrorMessage("Expected: enthalpy, Found: " + string(tag));
	for(j=1;j<=_n_variance;j++)
		for(i=1;i<=_n_csi;i++)
			fInput.read(reinterpret_cast < char * > (&_enthalpy_f[j][i]),sizeof(double));

	fInput.read(tag, SIZE);	// 8
	if (string(tag) != "mw")
		ErrorMessage("Expected: mw, Found: " + string(tag));
	for(j=1;j<=_n_variance;j++)
		for(i=1;i<=_n_csi;i++)
			fInput.read(reinterpret_cast < char * > (&_mw_f[j][i]),sizeof(double));

	fInput.read(tag, SIZE);	// 9
	if (string(tag) != "cp")
		ErrorMessage("Expected: cp, Found: " + string(tag));
	for(j=1;j<=_n_variance;j++)
		for(i=1;i<=_n_csi;i++)
			fInput.read(reinterpret_cast < char * > (&_cp_f[j][i]),sizeof(double));

	fInput.read(tag, SIZE);	// 10
	if (string(tag) != "lambda")
		ErrorMessage("Expected: lambda, Found: " + string(tag));
	for(j=1;j<=_n_variance;j++)
		for(i=1;i<=_n_csi;i++)
			fInput.read(reinterpret_cast < char * > (&_lambda_f[j][i]),sizeof(double));

	fInput.read(tag, SIZE);	// 11
	if (string(tag) != "mu")
		ErrorMessage("Expected: mu, Found: " + string(tag));
	for(j=1;j<=_n_variance;j++)
		for(i=1;i<=_n_csi;i++)
			fInput.read(reinterpret_cast < char * > (&_mu_f[j][i]),sizeof(double));

	fInput.read(tag, SIZE);	// 12
	if (string(tag) != "as")
		ErrorMessage("Expected: as, Found: " + string(tag));
	for(j=1;j<=_n_variance;j++)
		for(i=1;i<=_n_csi;i++)
			fInput.read(reinterpret_cast < char * > (&_as_f[j][i]),sizeof(double));
			
	// alpha
	for(j=1;j<=_n_variance;j++)
		for(i=1;i<=_n_csi;i++)
			_alpha_f[j][i] = _lambda_f[j][i]/_density_r[j][i]/_cp_f[j][i];				

	if(_number_of_species_to_extract > 0)
	{
		_mw_species_to_extract.resize(_number_of_species_to_extract+1);
		_w_f.resize(_number_of_species_to_extract+1);
		for(int k=1;k<=_number_of_species_to_extract;k++)
		{
			_mw_species_to_extract[k] = 0.;
			_w_f[k].resize(_n_variance+1);
			for(j=1;j<=_n_variance;j++)
				_w_f[k][j].resize(_n_csi+1);
		}
		
		int _count_species_found = 0;
		for(int k=1;k<=_number_of_species;k++)
		{
			double mw;
			double dummy;

			fInput.read(tag, SIZE);
			fInput.read(reinterpret_cast < char * > (&mw),sizeof(double));

			bool wFound = false;
			for(int w=1;w<=_number_of_species_to_extract;w++)
				if (_names_of_species_to_extract[w] == string(tag))
				{
					for(j=1;j<=_n_variance;j++)
						for(i=1;i<=_n_csi;i++)
							fInput.read(reinterpret_cast < char * > (&_w_f[w][j][i]),sizeof(double));
					
					_mw_species_to_extract[w] = mw;

					_count_species_found++;
					wFound = true;
					break;
				}

			if (_count_species_found == _number_of_species_to_extract)
				break;
			
			if (wFound == false)				// TODO
				for(j=1;j<=_n_variance;j++)
						for(i=1;i<=_n_csi;i++)
							fInput.read(reinterpret_cast < char * > (&dummy),sizeof(double));
		}
		
		for(int w=1;w<=_number_of_species_to_extract;w++)
			if (_mw_species_to_extract[w] == 0.)
				ErrorMessage("The following species is not available in the flamelet library: " + _names_of_species_to_extract[w]);
	}

	fInput.close();


	for(j=1;j<=_n_variance;j++)
		for(i=1;i<=_n_csi;i++)
		{
			if (_temperature_f_max < _temperature_f[j][i]) 
				_temperature_f_max = _temperature_f[j][i];

			if (_temperature_f_min > _temperature_f[j][i]) 
				_temperature_f_min = _temperature_f[j][i];

			if (_density_r_max < _density_r[j][i]) 
				_density_r_max = _density_r[j][i];

			if (_density_r_min > _density_r[j][i]) 
				_density_r_min = _density_r[j][i];
		}

	if ( _temperature_f_max <= 1.05*max(_temperature_f[1][1],_temperature_f[1][_n_csi]) )
		_cold = true;
}

void OpenSMOKE_PDF_Flamelet::GetMeanValues(const double csi, const double csiv2, vector<double>& extracted)
{
	int iCsi = 0;
	int jCsiv2 = 0;
	
	for(int j=2;j<=_n_variance;j++)
		if (csiv2 <= _variance_normal[j])
		{
			jCsiv2 = j-1;
			break;
		}

	for(int i=2;i<=_n_csi;i++)
		if (csi <= _csi[i])
		{
			iCsi = i-1;
			break;
		}
		
	double q   = (_csi[iCsi+1]-_csi[iCsi])*(_variance_normal[jCsiv2+1]-_variance_normal[jCsiv2]);
	double q11 = (_csi[iCsi+1]-csi) * (_variance_normal[jCsiv2+1]-csiv2);
	double q21 = (csi-_csi[iCsi])   * (_variance_normal[jCsiv2+1]-csiv2);
	double q12 = (_csi[iCsi+1]-csi) * (csiv2-_variance_normal[jCsiv2]);
	double q22 = (csi-_csi[iCsi])   * (csiv2-_variance_normal[jCsiv2]);
	
	// Temperature Favre [K]
	extracted[1]	= ( 	_temperature_f[jCsiv2][iCsi]*q11   + _temperature_f[jCsiv2][iCsi+1]*q21 +
	           				_temperature_f[jCsiv2+1][iCsi]*q12 + _temperature_f[jCsiv2+1][iCsi+1]*q22 ) / q;

	// Density Reynolds [kg/m3]
	extracted[2]	= (	_density_r[jCsiv2][iCsi]*q11   + _density_r[jCsiv2][iCsi+1]*q21 +
	          			_density_r[jCsiv2+1][iCsi]*q12 + _density_r[jCsiv2+1][iCsi+1]*q22 ) / q;
	          		
	// Absorption Coefficient Favre [1/m]
	extracted[3]	= (	_as_f[jCsiv2][iCsi]*q11   + _as_f[jCsiv2][iCsi+1]*q21 +
						_as_f[jCsiv2+1][iCsi]*q12 + _as_f[jCsiv2+1][iCsi+1]*q22 ) / q;         
						
	// Dynamic Viscosity [kg/m/s]
	extracted[4]	= (	_mu_f[jCsiv2][iCsi]*q11   + _mu_f[jCsiv2][iCsi+1]*q21 +
						_mu_f[jCsiv2+1][iCsi]*q12 + _mu_f[jCsiv2+1][iCsi+1]*q22 ) / q;   
						
	// Thermal diffusivity [m2/s]
	extracted[5]	= (	_alpha_f[jCsiv2][iCsi]*q11   + _alpha_f[jCsiv2][iCsi+1]*q21 +
						_alpha_f[jCsiv2+1][iCsi]*q12 + _alpha_f[jCsiv2+1][iCsi+1]*q22 ) / q;   	

	// Specific heat [J/kg/K]
	extracted[6]	= (	_cp_f[jCsiv2][iCsi]*q11   + _cp_f[jCsiv2][iCsi+1]*q21 +
						_cp_f[jCsiv2+1][iCsi]*q12 + _cp_f[jCsiv2+1][iCsi+1]*q22 ) / q;  																	  		
}

void OpenSMOKE_PDF_Flamelet::ExtractMeanValues(const double csi, const double csiv2, vector<double> &omegaFavre)
{
	int iCsi = 0;
	int jCsiv2 = 0;
	
	for(int j=2;j<=_n_variance;j++)
		if (csiv2 <= _variance_normal[j])
		{
			jCsiv2 = j-1;
			break;
		}

	for(int i=2;i<=_n_csi;i++)
		if (csi <= _csi[i])
		{
			iCsi = i-1;
			break;
		}
		
	double q   = (_csi[iCsi+1]-_csi[iCsi])*(_variance_normal[jCsiv2+1]-_variance_normal[jCsiv2]);
	double q11 = (_csi[iCsi+1]-csi) * (_variance_normal[jCsiv2+1]-csiv2);
	double q21 = (csi-_csi[iCsi])   * (_variance_normal[jCsiv2+1]-csiv2);
	double q12 = (_csi[iCsi+1]-csi) * (csiv2-_variance_normal[jCsiv2]);
	double q22 = (csi-_csi[iCsi])   * (csiv2-_variance_normal[jCsiv2]);
	
	for(int w=1;w<=_number_of_species_to_extract;w++)
		omegaFavre[w] = (	 _w_f[w][jCsiv2][iCsi]*q11		+ _w_f[w][jCsiv2][iCsi+1]*q21 +
	           				 _w_f[w][jCsiv2+1][iCsi]*q12	+ _w_f[w][jCsiv2+1][iCsi+1]*q22 ) / q;
}

void OpenSMOKE_PDF_Flamelet::Summary()
{
	cout << endl;
	cout << "OpenSMOKE_PDF_Flamelet" << endl;
	cout << "Index:                             " << _index << endl;
	cout << "File name:                         " << _name << endl;
	cout << "Mixture fraction points:           " << _n_csi << endl;
	cout << "Variance Mixture fraction points:  " << _n_variance << endl;
	cout << "Nominal strain rate:               " << _as << " Hz"  << endl;
	cout << "Stoic. scalar dissipation rate:    " << _chi_st << " Hz"  << endl;
	cout << "Temperature min:                   " << _temperature_f_min << " K" << endl;
	cout << "Temperature max:                   " << _temperature_f_max << " K" << endl;
	cout << "Density min:                       " << _density_r_min << " kg/m3" << endl;
	cout << "Density max:                       " << _density_r_max << " kg/m3" << endl;	
	if (_cold == true)	cout << "Cold flame" << endl;
	cout << endl;
}

void OpenSMOKE_PDF_Flamelet::ErrorMessage(const string error_message)
{
	cout << "OpenSMOKE_PDF_Flamelet Error" << endl;
	cout << "File name:     " << _name << endl;
	cout << "Error message: " << error_message << endl;
	exit(-1);
}

void OpenSMOKE_PDF_Flamelet::WarningMessage(const string warning_message)
{
	cout << "OpenSMOKE_PDF_Flamelet Warning" << endl;
	cout << "File name:     " << _name << endl;
	cout << "Warning message: " << warning_message << endl;
	cout << "Press enter to continue..." << endl;
	getchar();
}

	double OpenSMOKE_PDF_Flamelet::temperature_f_max()
	{
		return _temperature_f_max;
	};

	double OpenSMOKE_PDF_Flamelet::temperature_f_min()
	{
		return _temperature_f_min;
	};

	double OpenSMOKE_PDF_Flamelet::density_r_max()
	{
		return _density_r_max;
	};

	double OpenSMOKE_PDF_Flamelet::density_r_min()
	{
		return _density_r_min;
	};

	bool OpenSMOKE_PDF_Flamelet::cold()
	{
		return _cold;
	};
	
	double OpenSMOKE_PDF_Flamelet::density_r_fuel()
	{
		return _density_r[1][_n_csi];
	};

	double OpenSMOKE_PDF_Flamelet::density_r_oxidizer()
	{
		return _density_r[1][1];
	};
	
	double OpenSMOKE_PDF_Flamelet::enthalpy_f_fuel()
	{
		return _enthalpy_f[1][_n_csi];
	};

	double OpenSMOKE_PDF_Flamelet::enthalpy_f_oxidizer()
	{
		return _enthalpy_f[1][1];
	};
	
	double OpenSMOKE_PDF_Flamelet::temperature_f_fuel()
	{
		return _temperature_f[1][_n_csi];
	};
	
	double OpenSMOKE_PDF_Flamelet::temperature_f_oxidizer()
	{
		return  _temperature_f[1][1];
	};


	
