/*
 * HpgForms.cpp
 *
 *  Created on for cython 10 May 2021
 * by Andy
 */

#include <complex>
#include <cmath>
#include "HpgForms.h"


double pi=3.14159265358979323846264338327950288419716939937510;
dcmplx I2(0.0, 1.0);

dcmplx calculateHMatrixElement(int target, double Ip, dcmplx pz, double py,  dcmplx ts, dcmplx Eft, dcmplx Aft, double alpha){

	dcmplx Hpg;
	switch(target){
	case He :
		Hpg = HpgHe;
		break;
	case HeTheta :
		Hpg=HpgHeTheta;
		break;
	case Ne :
		Hpg = HpgNe;
		break;
	case Ar :
		Hpg = HpgAr;
		break;
	case ArEx_4S :
		Hpg = HpgArEx_4S;
		break;
	case Xe :
		Hpg = HpgXe;
		break;
	case N2 :
		Hpg = HpgN2Par;
		break;
	case N2Theta :
		Hpg = HpgN2Theta;
		break;
	case O2Theta :
		Hpg = HpgO2Theta;
		break;
	default :
		Hpg=1.;
		break;
	}
	return Hpg;
}