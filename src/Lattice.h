/*
 * Lattice.h
 *
 *  Created on: Mar 17, 2010
 *      Author: Ken
 */

#ifndef LATTICE_H_
#define LATTICE_H_
#include "Includes.h"
#include "Global.h"
class Lattice {
public:
	Lattice();
	virtual ~Lattice();
	double GetNNDist();
	double GetNN2Dist();
	int GetnNN();
	void LoadLatticePoints(string, string);
	void InitLatticeVect();
	void CalcRecipLattVects();
	void CalcRealLattVects();
	void AllocPoints();
	void CopyLattice(Lattice *);
	void GetCrossNNVectL(double *);
	double CalcVol(string);


	string name, symmetry;
	double a, b, c, alpha, beta, gamma;
	double vectors [3][3], recipVectors[3][3];
	int nPoints;
	double **points;
};

#endif /* LATTICE_H_ */
