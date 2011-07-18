/*
 * Timestep.h
 *
 *  Created on: Nov 18, 2010
 *      Author: ken
 */

#ifndef TIMESTEP_H_
#define TIMESTEP_H_
#include "Includes.h"
#include "Global.h"
#include "Group.h"
#include "Position.h"
#include "Velocity.h"
#include "SRList.h"


class Timestep {
public:
	Timestep();
	virtual ~Timestep();

	int step;
	Group all;
	int nUserGroups, nPosData, nVelData, nInPosActions;
	string *inPosAction;
	Group **userGroup;
	Group **posData, **velData, *firstPos, *firstVel, *curPos, *curVel;
	SRList sRList;
	Timestep *nextStep;

	void InputVel(string);
	void InputPos(string);
	void DoInputPosActions(Group *);
	void CompileAllVel();
	void CompileAllPos();
	bool GroupsAreEqual(Group*, Group*);
	void GetMaxMinAtomCoords();
	void BuildNearestNeighs(double);
	void BuildNeighsInRange(double, double);
	Group * FindGroup(string);
	Group * FindPos(string);
	void CreatePosArray();
	void CreatePosInputAction(string);
	void DeletePosData(string);
	void CreateVelArray();
	void ReadVelInLammpsVelFile(Group *, Velocity *);
	void ReadPosInLammpstrjFile(Group *, Position *);

	//void CreateUserGroup();
	void Clean();
	void CalcLatticeParameter(Lattice *);

	// User group creation (in TS_Groups.cpp)
	void AddUserGroup(string, string);
	void GenGroups();
	void GroupFromType(Group *, Group *, string);
	void GroupFromNN(Group *, Group *, string);
};

#endif /* TIMESTEP_H_ */
