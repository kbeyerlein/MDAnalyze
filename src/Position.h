/*
 * Position.h
 *
 *  Created on: Jan 3, 2010
 *      Author: Ken
 */

#ifndef POSITION_H_
#define POSITION_H_

#include "Includes.h"
#include "Global.h"
#include "Lattice.h"
#include <gsl/gsl_math.h>

class Position
{
public:
	Position();
	~Position(void);

	int dataPos;
	string file, fileType;
	gsl_vector *center;
	double dim[3][2], sideL[3];
	Distrib *rad, *rdf;
	Lattice lattice;
	char boundary[3];


	void GetDimensions();
	double GetMaxD();
	//void CalcAvgLatticeParam();
	void CalcLatticeVectors();
	void FindApproxLattPlanes(gsl_vector **, int , double *, Lattice *);
	bool SameBoxInfo(Position *);




};

#endif /* POSITION_H_ */
