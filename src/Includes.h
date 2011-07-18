/*
 * Includes.h
 *
 *  Created on: Aug 24, 2009
 *      Author: Ken
 */

#ifndef INCLUDES_H_
#define INCLUDES_H_
#include <iostream>
#include <iomanip>
//#include <stdlib.h>
#include <fstream>
#include <sstream>
#define _USE_MATH_DEFINES
#include <cmath>

#include <gsl/gsl_blas.h>

#define getmax(a,b) a>b?a:b
#define getmin(a,b) a<b?a:b
#define PI M_PI
#define PREC .00001f
#define SRPREC .001f
#define DPREC .001f
#define VPREC .001f
using namespace std;

struct Param{
	string name;
	double val;
};
/*
struct Lattice{
	string name, symmetry;
	double a, b, c, alpha, beta, gamma;
	double vectors [3][3], recipVectors[3][3];
	int nPoints;
	double **points;
};
*/
struct Atom{
	int index;
	int type;
	double r, v;
	gsl_vector *pos, *vel;
	int n_NN;
	Atom **NN;
	// double mass, charge;
};
struct Basis{
	string name;
	int nAtoms;
	Atom *atoms;
};
struct Distrib{
	string name, xName, xUnits, yName, yUnits;
	int n;
	double step,minX,maxX;
	double *x, *y;
};
struct Range{
	string name;
	double x0,xf,dx;
};


#endif /* INCLUDES_H_ */
