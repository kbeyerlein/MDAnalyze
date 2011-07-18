/*
 * Displace.h
 *
 *  Created on: Jan 11, 2010
 *      Author: Ken
 */

#ifndef DISPLACE_H_
#define DISPLACE_H_

#include "Includes.h"
#include "Global.h"
#include "Group.h"
class Displace {
public:
	Displace();
	virtual ~Displace();
	void Allocate();
	void CalcDisplacements(double, double);
	double GetMaxDisplaceSq();
	void CalcPMSD();
	void CalcAvgMSD();
	void CalcRadMSD(gsl_vector *);
	void CalcRadDisp(gsl_vector *);
	void OutputMSD(string, string );
	void Zero();

	Group *pos_i, *pos_f;
	gsl_vector **u;
	double *uSq, msd;
	int nAtoms;
	string specifyMode, groupName;
	int initial, final;
	Distrib *pMSD, *radMSD, *radDisp;
	bool scaled;
};

#endif /* DISPLACE_H_ */
