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
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

#ifndef OpenSMOKE_ChiDistribution_H
#define OpenSMOKE_ChiDistribution_H

class OpenSMOKE_ChiDistribution
{

public:

	OpenSMOKE_ChiDistribution();	// Constructor
	void BuildGrid();
	void AssignScalarDissipationRates(const std::vector<double> Xst);
	void AssignMeanScalarDissipationRate(const double XstMean);
	
	double GetMeanValue(const std::vector<double> &values);
	std::vector<double> ExtractMeanValue(const std::vector < std::vector<double> > &values);

public:

	void SetSigma(const double sigma);
	void SetLowerLimit(const double user_Xl_min);
	void SetHigherLimit(const double user_Xl_max);
	void SetNumberOfPoints(const int nIntervals);
	void SetFixedPointRatio();
	void SetAccurateCalculation();
	void SetWeightThreshold(const double weight_threshold);

	inline int iXinf() { return _iXinf;}  	// Werte zurückgeben
	inline int iXsup() { return _iXsup;}	// Werte zurückgeben

private:

	std::string _name;								// library name
	
	double _sigmaSquared;
	double _teta_coeff;
	double _mu;
	double _Xq;									// extinction scalar dissipation rate

private:

	int		_nIntervals;
	bool	_iFixedRatioStep;
	double	_sigma;
	double	_user_Xl_min;
	double	_user_Xl_max;
	double	_weight_threshold;
	bool	_iAccurate;
	
private:

	int _iXmax;
	int _iXinf;		// Funktion für Rückgabe verfügbar
	int _iXsup;		// Funktion für Rückgabe verfügbar
	int _index_inf;
	int _index_sup;

	std::vector<double>	_Xst;
	std::vector<double>	_Xl;
	std::vector<double>	_Xc;
	std::vector<double>	_teta;
	std::vector<double>	_erfteta;
	std::vector<double>	_weights;
	std::vector<double>	_ratioX;
	std::vector<double>	_y;
	std::vector<int>		_iX;

private:

	void PrepareInterpolation();
	
	void Summary();
	void SummaryDebug(const double XstMean);
	void SummaryDebug();

	void ErrorMessage(const std::string error_message);
	void WarningMessage(const std::string warning_message);
};

#endif // OpenSMOKE_PDF_Flamelet_Library_H
