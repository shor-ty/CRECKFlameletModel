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
#include <iomanip>
#include <math.h>

using namespace std;

#include "OpenSMOKE_ChiDistribution.hpp"
#include "OpenSMOKE_External_Functions.hpp"

OpenSMOKE_ChiDistribution::OpenSMOKE_ChiDistribution()
{
	_name = "[not assigned]";
	
	_sigma			= 1.31;
	_user_Xl_min		= 1.e-6;
	_user_Xl_max		= 5.e2;
	_nIntervals		= 40;
	_iFixedRatioStep	= true;
	_weight_threshold	= 1.e-10;
	_iAccurate		= false;

	_sigmaSquared		= _sigma*_sigma;
	_teta_coeff		= 1./sqrt(2.)/_sigma;
}

void OpenSMOKE_ChiDistribution::BuildGrid()
{
	double _alfa;
	if (_iFixedRatioStep == true)
	{
		double deltax;
		_alfa = 1.01;
		double ratio = 2.;
		int outer_iterations = 100;
		int inner_iterations =  30;

		for(int k=1;k<=outer_iterations;k++)
		{
			double sum = 1.;
			for(int j=1;j<=_nIntervals-3;j++)
				sum += pow(_alfa,j); 
			deltax = 2.*_user_Xl_min*sum;
			
			if (deltax < (_user_Xl_max-_user_Xl_min))
				_alfa *= ratio;
			else
			{
				double alfa_a = _alfa/ratio;
				double alfa_b = _alfa;
				
				for (int n=1;n<=inner_iterations;n++)
				{
					_alfa = 0.50*(alfa_a+alfa_b);
					
					double sum = 1.;
					for(int j=1;j<=_nIntervals-3;j++)
						sum += pow(_alfa,j); 
					deltax = 2.*_user_Xl_min*sum;

					if (deltax < (_user_Xl_max-_user_Xl_min))
						alfa_a = _alfa;
					else
						alfa_b = _alfa;

					if (fabs(deltax-(_user_Xl_max-_user_Xl_min))/deltax < 0.0001)
					{	
						_alfa = 0.50*(alfa_a+alfa_b);
						break;
					}
				}
				break;
			}
		}

		_Xl.resize(_nIntervals+1);
		_Xl[1] = 0.;
		_Xl[2] = _user_Xl_min;
		for(int j=3;j<=_nIntervals;j++)
			_Xl[j] = _Xl[j-1] + 2.*_user_Xl_min*pow(_alfa, j-3.);
	}

	else
	{
		_alfa = pow(_user_Xl_max/_user_Xl_min, 1./double(_nIntervals-2.));
		_Xl.resize(_nIntervals+1);
		_Xl[1] = 0.;
		_Xl[2] = _user_Xl_min;
		for(int j=3;j<=_nIntervals;j++)
			_Xl[j] = _Xl[j-1]*_alfa;
	}

	double difference = fabs(_Xl[_nIntervals]-_user_Xl_max)/_user_Xl_max*100.;
	if (difference > 0.1)
	{
	ErrorMessage("Please change the Chi Distribution input parameters...");

	cout << "-----------------------------------------------------------" << endl;
	cout << "               OpenSMOKE_ChiDistribution               " << endl;
	cout << "-----------------------------------------------------------" << endl;
	cout << " Increment ratio:      " << _alfa << endl;
	cout << " User defined Xl_min:  " << _user_Xl_min		<< " [Hz]" << endl;
	cout << " User defined Xl_max:  " << _user_Xl_max		<< " [Hz]" << endl;
	cout << " Calculated Xl_max:    " << _Xl[_nIntervals]	<< " [Hz]" << endl;
	cout << " Difference (%):       " << difference					   << endl;
	cout << endl;
	}
	else
	cout << "hier";

	// Memory allocation
	_Xc.resize(_nIntervals+1);
	_teta.resize(_nIntervals+1);
	_erfteta.resize(_nIntervals+1);
	_weights.resize(_nIntervals+1);
	_y.resize(_nIntervals+1);

	// Build centered grid
	for(int i=1;i<=_nIntervals-1;i++)
		_Xc[i] = 0.50*(_Xl[i+1]+_Xl[i]);
	_Xc[_nIntervals] = _Xl[_nIntervals]*2.;

	_teta[1]	= -1.e16;	// -Inf
	_erfteta[1]	= -1.;		//

}

void OpenSMOKE_ChiDistribution::SetSigma(const double sigma)
{
	if (sigma < 0.50 || sigma > 5)
		ErrorMessage("Sigma must be between 0.50 and 5 (default: 1.31)");
	_sigma			= sigma;
	_sigmaSquared		= _sigma*_sigma;
	_teta_coeff		= 1./sqrt(2.)/_sigma;
}

void OpenSMOKE_ChiDistribution::SetLowerLimit(const double user_Xl_min)
{
	if (user_Xl_min < 1.e-10 || user_Xl_min > 1.e-1)
		ErrorMessage("Lower limit of scalar dissipation grid must be between 1.e-10 and 1.e-1 (default: 1.e-6)");
	_user_Xl_min = user_Xl_min;
}

void OpenSMOKE_ChiDistribution::SetHigherLimit(const double user_Xl_max)
{
	if (user_Xl_max < 1. || user_Xl_max > 1.e6)
		ErrorMessage("Lower limit of scalar dissipation grid must be between 1.e0 and 1.e6 (default: 1.e3)");
	_user_Xl_max = user_Xl_max;
}

void OpenSMOKE_ChiDistribution::SetNumberOfPoints(const int nIntervals)
{
	if (nIntervals < 10 || nIntervals > 60)
		ErrorMessage("Number of grid points must be between 10 and 60 (default: 40)");
	_nIntervals = nIntervals;
}

void OpenSMOKE_ChiDistribution::SetFixedPointRatio()
{
	_iFixedRatioStep = false;
}

void OpenSMOKE_ChiDistribution::SetAccurateCalculation()
{
	_iAccurate = true;
}
	
void OpenSMOKE_ChiDistribution::SetWeightThreshold(const double weight_threshold)
{
	if (weight_threshold < 1.e-16 || weight_threshold > 1.e-4)
		ErrorMessage("Weight threshold must be between 1.e-16 and 1.e-4 (default: 1.e-10)");
	_weight_threshold = weight_threshold;
}
	
void OpenSMOKE_ChiDistribution::AssignScalarDissipationRates(const vector<double> Xst)
{
	_Xst = Xst;						
	_Xq = _Xst[_Xst.size()-1];		

	if (_Xc[_nIntervals] <= _Xq)
		ErrorMessage("The list of scalar dissipation rates is too short to cover the extinction limits...");

	for(int i=1;i<=_nIntervals;i++)
		if (_Xc[i]>_Xq)
		{
			_iXmax			= i-1;
			break;
		}
	_iX.resize(_iXmax+1);
	_ratioX.resize(_iXmax+1);

	_index_inf = 1;
	_index_sup = _nIntervals;
	_iXinf = 1;
	_iXsup = _Xst.size()-1;

	PrepareInterpolation();

	Summary();
}

void OpenSMOKE_ChiDistribution::PrepareInterpolation()
{
	for(int i=1;i<=_iXmax;i++)
		for(int j=2;j<=int(_Xst.size()-1);j++)
		{
			if (_Xc[i] <= _Xst[j])
			{	
				_iX[i] = j-1;
				_ratioX[i] = (_Xc[i]-_Xst[j-1])/(_Xst[j]-_Xst[j-1]);
				break;
			}
		}
}

void OpenSMOKE_ChiDistribution::AssignMeanScalarDissipationRate(const double XstMean)
{
	_mu   = log(XstMean)-0.50*_sigmaSquared;

	// _teta[1] = -Inf 
	for(int i=2;i<=_nIntervals;i++)
		_teta[i] = _teta_coeff*(log(_Xl[i]/XstMean)+0.50*_sigmaSquared);

	// _erfteta[1] = -1
	for(int i=2;i<=_nIntervals;i++)
		_erfteta[i] = sun_erf(_teta[i]);

	// weights
	for(int i=1;i<=_nIntervals-1;i++)
		_weights[i] = 0.50*(_erfteta[i+1]-_erfteta[i]);
	_weights[_nIntervals] = 0.50*(1.-_erfteta[_nIntervals]);

	if (_iAccurate == false)
	{
		_index_inf = 1;
		for(int i=1;i<=_nIntervals;i++)
		{
			if (_weights[i]>_weight_threshold)	
			{
				_index_inf = i;
				break;
			}
		}

		_index_sup = _nIntervals;
		for(int i=_index_inf+1;i<=_nIntervals;i++)
		{
			if (_weights[i]<=_weight_threshold)	
			{
				_index_sup = i-1;
				break;
			}
		}

		_iXinf = _iX[_index_inf];
		_iXsup =_Xst.size()-1;
		if (_index_sup <= _iXmax)
			_iXsup = _iX[_index_sup]+1;
	}

//	SummaryDebug(XstMean);
}



double OpenSMOKE_ChiDistribution::GetMeanValue(const vector<double> &values)
{
	for(int i=1;i<=_iXmax;i++)	
		_y[i] = values[_iX[i]]+_ratioX[i]*(values[_iX[i]+1]-values[_iX[i]]);
	for(int i=_iXmax+1;i<=_nIntervals;i++)	
		_y[i] = values[values.size()-1];

	double sum=0.;
	for(int i=_index_inf;i<=_index_sup;i++)
		sum += _y[i]*_weights[i];

//	SummaryDebug();
	
	return sum;
}

vector<double> OpenSMOKE_ChiDistribution::ExtractMeanValue(const vector < vector<double> > &values)
{
	vector<double> sum;
	sum.resize(values[1].size());

	for(int w=1;w<=int(values[1].size()-1);w++)
	{
		for(int i=1;i<=_iXmax;i++)	
			_y[i] = values[_iX[i]][w]+_ratioX[i]*(values[_iX[i]+1][w]-values[_iX[i]][w]);
		for(int i=_iXmax+1;i<=_nIntervals;i++)	
			_y[i] = values[values.size()-1][w];

		sum[w] = 0.;
		for(int i=_index_inf;i<=_index_sup;i++)
			sum[w] += _y[i]*_weights[i];
	}

	return sum;
}

void OpenSMOKE_ChiDistribution::ErrorMessage(const string error_message)
{
	cout << "OpenSMOKE_ChiDistribution Error" << endl;
	cout << "Error message: " << error_message << endl;
	exit(-1);
}

void OpenSMOKE_ChiDistribution::WarningMessage(const string warning_message)
{
	cout << "OpenSMOKE_ChiDistribution Warning" << endl;
	cout << "Warning message: " << warning_message << endl;
	cout << "Press enter to continue..." << endl;
	getchar();
}

void OpenSMOKE_ChiDistribution::Summary()
{
	cout << endl;
	cout << " Number of intervals:                     " << _nIntervals			<< endl;
	cout << " Sigma:                                   " << _sigma				<< endl;
	cout << " Extinction scalar dissipation rate (Xq): " << _Xq					<< " Hz" << endl;
	cout << " Maximum scalar dissipation rate (Xc):    " << _Xc[_nIntervals]	<< " Hz" << endl;
	cout << endl;
	cout << "-----------------------------------------------------------" << endl;
	cout << " Point            Xl               Xc        Index   " << endl;
	cout << "-----------------------------------------------------------" << endl;
	for(int i=1;i<=_iXmax;i++)
		cout << setw(4) << i << setw(16) << _Xl[i] << setw(16) << _Xc[i] << setw(8) << _iX[i] << endl;
	for(int i=_iXmax+1;i<=_nIntervals;i++)
		cout << setw(4) << i << setw(16) << _Xl[i] << setw(16) << _Xc[i] << setw(8) << _Xst.size()-1 << endl;
	cout << endl;
}

void OpenSMOKE_ChiDistribution::SummaryDebug()
{
	cout << endl;
	cout << "OpenSMOKE_ChiDistribution" << endl;
	cout << "-----------------------------------------------------------" << endl;
	cout << "Point   Xc              W                                  " << endl;
	cout << "-----------------------------------------------------------" << endl;
	for(int i=1;i<=_nIntervals;i++)
		cout << i << "\t" << _Xc[i] << "\t" << _weights[i] << "\t" << _y[i]*_weights[i] << endl;

	cout << endl;
}

void OpenSMOKE_ChiDistribution::SummaryDebug(const double XstMean)
{
	cout << endl;
	cout << "OpenSMOKE_ChiDistribution" << endl;
	cout << "-----------------------------------------------------------" << endl;
	cout << "Point   Xc               W                                 " << endl;
	cout << "-----------------------------------------------------------" << endl;
	for(int i=1;i<=_nIntervals;i++)
		cout << i << "\t" << _Xc[i] << "\t" << _weights[i] << endl;
	cout << endl;
}
