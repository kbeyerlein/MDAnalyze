/*
 * Global.h
 *
 *  Created on: Mar 15, 2010
 *      Author: Ken
 */

#ifndef GLOBAL_H_
#define GLOBAL_H_
#include "Includes.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>

double Dot(double *, double *);
void Cross(double *, double *, double *);
void Cross(gsl_vector*, gsl_vector*, gsl_vector*);
double VectorL(gsl_vector*);
void InitBasis(Basis *);
void Normalize(double *);
void ScaleVect(double *, double);
void CalcCenterOfMass_gsl(gsl_vector**, int,double *, gsl_vector*);
void OuterProduct_gsl(gsl_vector *, gsl_vector *, gsl_matrix *);
double CalcJacobian_gsl(gsl_matrix *);
void CalcDistPBC_gsl(gsl_vector *, gsl_vector *, gsl_vector *, char *, double *, double *);

void QuickSortImproved(double *, double *, int , int);
void SortByInsertion(double *, double *, int , int );
int Partition(double *, double *, int , int );
void CopyAtom(Atom *, Atom *);
void ZeroAtom(Atom *);
void AllocateArray(double**, int);
void InitDistrib(Distrib *);
void AllocAtomNN(Atom *);
gsl_vector * ArrayColToGSLVector(double **, int, int);
void OutputDistrib(Distrib*, string, string, bool, bool, string);
void OutputMultDistribs(Distrib **dist, int, string, string, bool, bool, bool, string);
void AllocDistrib(Distrib **, string, string, string, double, double, double);
void AllocDistrib(Distrib *, string, string, string, double, double, double);
void ScaleDistrib(Distrib *, string, double);
void CleanDistrib(Distrib *);
void FindBestPlaneOfPoints(gsl_vector**, double *, int, gsl_vector *, double *);
void OutputStrainTensor(gsl_matrix **, int, string, string, string);

// Global variables and Functions for Ranges
//#ifndef SRANGE
//#define SRANGE
//int nSpaceRange=0;
//Range **sRange=0;
//#endif


#endif /* GLOBAL_H_ */
