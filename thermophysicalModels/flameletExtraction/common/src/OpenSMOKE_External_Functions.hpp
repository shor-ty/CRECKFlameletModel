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
#include <math.h>

using namespace std;

#ifndef OpenSMOKE_External_Functions_H
#define OpenSMOKE_External_Functions_H

double sun_erf(double x);
double sun_erfc(double x);

class FractionalExponents
{

public:

	FractionalExponents(const double CStar_, const double delta_, const double alfa_);
	
	double gamma(const double &c, const double lambda);

private:

	double alfa;
	double H;
	double K;
	double CStar;
	double delta;
};


enum chiPDF_kinds		{	CHI_PDF_DIRAC, CHI_PDF_LOG_NORMAL};

enum nucleation_models 		{	NUCLEATION_NONE, NUCLEATION_LIU_2001, NUCLEATION_LIU_2002, NUCLEATION_MOSS_1999,
					NUCLEATION_WEN_2003, NUCLEATION_LINDSTEDT_1994, NUCLEATION_LEUNG_1991};

enum growth_models 		{	GROWTH_NONE, GROWTH_LIU_2001, GROWTH_LIU_2002, GROWTH_MOSS_1999,
					GROWTH_WEN_2003, GROWTH_LINDSTEDT_1994, GROWTH_LEUNG_1991};

enum aggregation_models 	{	AGGREGATION_NONE, AGGREGATION_SMOLUCHOWSKI, AGGREGATION_MOSS};
												
enum oxidation_models 		{	OXIDATION_NONE, OXIDATION_LEE, OXIDATION_NSC, OXIDATION_NEOH};

enum surface_functions		{	SURFACE_FUNCTION_LINEAR, SURFACE_FUNCTION_SQUARE_ROOT };


class OpenSMOKE
{
public:
	// Constants
	static const double pi;
	static const double Nav_kmol ;
	static const double rho_soot;
	static const double kB;   // Boltzmann constant [J/K]


	// Element molecular weights
	static const double mw_c;
	static const double mw_h;
	static const double mw_o;
	static const double mw_n;

	// Species Molecular weights
	static const double mw_c2h2;
	static const double mw_o2;
	static const double mw_oh;
	static const double mw_co;
	static const double mw_h2;

	// Constants
	static const double _36pi_to_1_over_3;
	static const double _8kB_over_pi;
	static const double _8_over_pi;
	static const double _4pi;
	static const double _6pi;
	static const double _16pi;
	static const double _sqrt2_over_2;
	static const double _4_over_sqrt2 ;
};

namespace Liu_2001
{
	const	double n_nuclei_		= 90000.;
	const	double A_nucleation_	= 30.0;
	const	double T_nucleation_	= 20643.;
	const	double A_growth_		= 12000.;
	const	double T_growth_		= 12083.;
	const	surface_functions	surface_function_ = SURFACE_FUNCTION_LINEAR;
}

namespace Liu_2002
{
	const	double n_nuclei_		= 700.;
	const	double A_nucleation_	= 2.85;
	const	double T_nucleation_	= 16103.;
	const	double A_growth_		= 42000.;
	const	double T_growth_		= 10064.;
	const	surface_functions	surface_function_ = SURFACE_FUNCTION_SQUARE_ROOT;
}	

namespace Moss_1999
{
	const	double n_nuclei_		= 12.;
	const	double A_nucleation_	= 54.0;
	const	double T_nucleation_	= 21100.;
	const	double A_growth_		= 11700.;
	const	double T_growth_		= 12100.;
	const	surface_functions	surface_function_ = SURFACE_FUNCTION_LINEAR;
}		

namespace Wen_2003
{
	const	double n_nuclei_		= 100.;
	const	double A_nucleation_	= 54.0;
	const	double T_nucleation_	= 21100.;
	const	double A_growth_		= 9000.6;
	const	double T_growth_		= 12100.;
	const	surface_functions	surface_function_ = SURFACE_FUNCTION_LINEAR;
}

namespace Lindstedt_1994
{
	const	double n_nuclei_		= 60.;
	const	double A_nucleation_	= 210.;
	const	double T_nucleation_	= 21100.;
	const	double A_growth_		= 18000.;
	const	double T_growth_		= 12100.;
	const	surface_functions	surface_function_ = SURFACE_FUNCTION_LINEAR;
}

namespace Leung_1991
{
	const	double n_nuclei_		= 32.;
	const	double A_nucleation_	= 635.0;
	const	double T_nucleation_	= 21100.;
	const	double A_growth_		= 144000.;
	const	double T_growth_		= 12100.;
	const	surface_functions	surface_function_ = SURFACE_FUNCTION_SQUARE_ROOT;
}		

namespace Lee
{
	const double A	= 8903.51;
	const double T	= 19778.;
}

namespace Neoh
{
	const double eta	= 0.04;
	const double A		= 105.81;
}

namespace NSC
{
	const double A	= 20.;
	const double TA	= 15098.;
	const double B	= 4.46e-3;
	const double TB	= 7650.;
	const double C	= 1.510e5;
	const double TC	= 48817.;
	const double D	= 21.3;
	const double TD	= 2063.;
}

namespace Smoluchowski
{
	const double gamma	= 7.94e-40;
}

namespace Moss
{
	const double gamma = 1.e9/(OpenSMOKE::Nav_kmol*OpenSMOKE::Nav_kmol);
}

double free_path(const double mu, const double P_atm, const double T, const double MW);


#endif // OpenSMOKE_External_Functions_H
	
