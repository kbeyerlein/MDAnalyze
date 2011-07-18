/*
 * TimeSeries.h
 *
 *  Created on: Jan 11, 2010
 *      Author: Ken
 */

#ifndef TIMESERIES_H_
#define TIMESERIES_H_
#include "Global.h"
#include "Displace.h"
#include "Timestep.h"
#include "SRList.h"
#include <gsl/gsl_linalg.h>

struct AtomStrain{
	int nA, t0, tf;
	Group *g_i;
	gsl_matrix **s;
	int *index;
};


class TimeSeries {
public:
	TimeSeries();
	virtual ~TimeSeries();

	//TODO Unsure functions
	void CalcInstantLattice(int);

	//TODO Unlikely class variables;
	Lattice *lattice, *relaxedLatt;
	//TODO Unsure class variables
	Displace *avgDynDisp;

	//gsl_matrix **E, **e;
	AtomStrain *F;
	int nStrTen;
	//Class variables
	string name;
	int nTimeSteps, nDynDisp;
	Timestep *firstStep, **step;
	Displace *staticDisp, **dynDisp, *userDisp;
	//TODO Consider rewriting avg to be a timestep so it can also have groups
	Group *avg;
	Distrib *velCorrFn, *msd_t;

	int nTimeRange;
	Range **tRange;

	SRList sRList;

	//File parsing functions
	void FindPosInLammpstrjFile(string, string, char *);
	void FindVelInLammpsVelFile(string, string);
	void InsertStepInSeq(Timestep **, Timestep **, int);
	void CreateSeriesArray();
	//Gets from array
	int GetSmallConstTimeStep(int, int ,string);
	int GetNumPosInRange(int, int);
	int GetNumPosInRange(int, int, int);
	int GetNumVelInRange(Range *);
	int GetMaxTimestep();
	Timestep *FindTimestep(int);
	//Velocity correlation functions
	void CalcVelCorrFunc( string, string);
	double CalcAvgVelAutoCorr(string, Timestep *, Timestep *);
	//Avg Position
	void AllocAvg(string);
	void CalcAvgPos(string, int, int);
	//Atom Displacements
	void CalcDisplace(Displace **, string, string , bool, int , int);
	void AssignDispPositions(Displace *);
	void CalcStaticDisp(string , bool, int , int, int);
	void CalcDynamicDisp(string, bool, int, int, int);
	void CalcMSD_t(string, int, int, int);
	//Atom Local Deformation Strain Tensor
	void CalcLocalDefGradTensor(gsl_vector **, gsl_vector **, gsl_vector *, int, gsl_matrix *);
	void CalcStrainTensor(Group *, Group *, double, string, AtomStrain *);
	void CalcStrainTensor(Timestep *, Timestep *,string,  double, string, gsl_matrix *);
	void ApplyHeavysideStepWeights(gsl_vector *, gsl_vector **, int, double);
	void CalcCGStrainTensor(gsl_matrix *, gsl_matrix *, string);
	void DecompDefGradTensor(gsl_matrix *, gsl_matrix*, gsl_matrix*, gsl_matrix*);
	void CalcAtomicStrainTensors(string, int, int, int, double , string );
	void CalcCumulativeStrainTensor(string, int , int, int, double, string);
	void OutputEigenVectors(string, string, string, gsl_vector_complex*, gsl_matrix_complex *);
	void OutputEigenVectors(string, string, string, gsl_vector*, gsl_matrix *);
	void CalcLGStrainTensor(int , AtomStrain *);
	void CalcEAStrainTensor(int , AtomStrain *);
	void CalcAvgAtomStrain(int, gsl_matrix *);
	void CleanAtomTensor(AtomStrain *);
	void OutputStrainTensor(AtomStrain *, string , string , string);
	// Special RDF Calcs
	void CalcAvgAtomDistInRange(string, int, double, double, double);
	void CalcNN_t(Distrib *, string, string, string, string);
	// TimeRange Creation and finding
	Range* CreateNewTimeRange();
	Range* FindTimeRange(string );

	// Creation/Annihilation of Timestep from data loaded in other Timesteps
	void MovePosToNewTime(int, string, int, string);
	// Creating user defined group
	void DefineUserGroup(string, string, string);
	// Analyze the radial positions
	void CalcRadius(string, string, string, string);
};

#endif /* TIMESERIES_H_ */
