/*
 * Group.h
 *
 *  Created on: Nov 18, 2010
 *      Author: ken
 */

#ifndef GROUP_H_
#define GROUP_H_

#include "Includes.h"
#include "Global.h"
#include "Position.h"
#include "Velocity.h"
#include "SRList.h"
class Group {
public:
	Group();
	virtual ~Group();

	string name, type;
	int nAtoms;
	Atom **atom;
	Position *p;
	Velocity *v;
	Group *next;
	double neighR0, neighRc;
	Group *neighG;
	int nSubDomX, nSubDomY, nSubDomZ;
	Group ***subDomain;
	string genCommand;

	void SortByIndex();
	void AllocAtoms();
	void AllocPos();
	void AllocVel();
	void AllocAtomPos();
	void AllocAtomVel();
	void CopyAtomInfo(Group *);
	void CopyAtomPos(Group *);
	void CopyAtomVel(Group *);
	void CalcVelDistrib();
	void CalcVelMags();
	void CalcPosMags(gsl_vector *);
	double GetMaxVel();
	double GetMaxR();
	double GetMinR();
	double GetAvgR();
	double GetVarR(double);
	double CalcRadiusOfGyration(gsl_vector *);
	void ScalePos(gsl_vector *);
	void TransPos(gsl_vector *);
	void RotatePos(gsl_matrix *, gsl_vector *);
	void RemapPos(int, double, double);
	void CalcCenter();
	int GetTotNumNN();
	void AllocNN();
	void CreateRadialDistrib(gsl_vector*, double);
	void OutputNumNNLammpstrjFile(string , string , int);
	void OutputNNVectorLammpstrjFile(string, string, int);
	void CalcRDF(double);
	void CalcRDF(Distrib *);
	void CalcDensityVsR(Distrib *);
	void OutputVolChangeMap(double *, string, string, int ,int);
	Atom* FindAtom(int );
	void OutputStrainTensor(gsl_matrix **, int, string, string, string);
	void BuildNeighsInRange(Group*, double, double);
	void CalcAvgAtomDistInRange(Group *, gsl_vector*, double, double, double);
	void CalcAvgNeighDistVsR(Distrib *);
	void CalcNumNNDist(Distrib *, double);
	void AllocSubDomains(int, int, int);
	void CreateSubDomains(int*);
	void BuildNeighsInRange_fast(Group *, double, double);
	void OutputPosition(string, string, string);
	bool CheckAtomsInBox();
	// Group maintenance and cleaning
	void ShrinkGroup();
	//Error Messages
	void NoGroupBoolean(string);
};

#endif /* GROUP_H_ */
