/*
 * TS_Groups.cpp
 *
 *  Created on: Jun 16, 2011
 *      Author: ken
 */

#include "Timestep.h"

void Timestep::GenGroups()
{
	CompileAllPos();
	CompileAllVel();
	for (int i=0;i<nUserGroups;i++){
		// The commands used to generate a user group can be divided into 3 subcommands:
		// 		(1) param: this defines the parameter which is used to define the group.
		// 			(i.e. index, atomType, position, velocity ...)
		//		(2) boolean: this is a word which defines whether the specified range is inclusive, exclusive or restrictively equal.
		//			(i.e. inVol, outVol, equal ...)
		// 		(3) def: this is the defining values used to find the atoms in the group.
		//			(i.e. volName, typeRange, indexRange ...)
		stringstream comm(userGroup[i]->genCommand);
		string param, boolean, def;
		comm >> param;
		if (param=="type"){
			getline(comm,def);
			GroupFromType(&all, userGroup[i], def);
		}
		else if (param=="nNN"){
			getline(comm,def);
			GroupFromNN(&all, userGroup[i], def);
		}
		else {
			cout<<"Error creating userGroup: "<<userGroup[i]->name<<", the parameter type "<< param<<" is not supported.\n";
			exit(0);
		}
	}
}
// Note: Do not destroy or generate new user groups themselves only compile the list of Atom*
//       (this will avoid problems with multiple InputPos commands during the same run)
void Timestep::AddUserGroup(string n, string comm)
{
	// Save pointers and reallocate userGroup
	Group **tempUGroups=new Group* [nUserGroups+1];
	for (int i=0;i<nUserGroups;i++){
		tempUGroups[i]=userGroup[i];
	}
	if (userGroup){
		delete [] userGroup;
	}
	nUserGroups++;
	userGroup=new Group * [nUserGroups];
	for (int i=0;i<nUserGroups-1;i++){
		userGroup[i]=tempUGroups[i];
	}

	//Create new group and give it the name and command
	userGroup[nUserGroups-1]=new Group();
	userGroup[nUserGroups-1]->name=n;
	userGroup[nUserGroups-1]->genCommand=comm;
	userGroup[nUserGroups-1]->type="user";

	delete [] tempUGroups;
}
void Timestep::GroupFromType(Group *parent, Group *child, string genComm)
{
	stringstream comm(genComm);
	string boolean;
	comm>>boolean;
	int groupType;
	if (boolean=="equal"){
		comm>>groupType;
	}
	else{
		child->NoGroupBoolean(boolean);
	}
	child->nAtoms=parent->nAtoms;
	child->AllocAtoms();
	int nCAtoms=0;
	for (int i=0; i<parent->nAtoms;i++){
		if (boolean=="equal"){
			if (parent->atom[i]->type==groupType){
				child->atom[nCAtoms]=parent->atom[i];
				nCAtoms++;
			}
		}
		else{
			child->NoGroupBoolean(boolean);
		}
	}
	child->ShrinkGroup();
	child->p=parent->p;
	child->v=parent->v;
	cout<<"Created group with "<<child->nAtoms<<endl;
}
void Timestep::GroupFromNN(Group *parent, Group *child, string genComm)
{
	//Command for defining group consists of: (boolean) (range) (spaceRangeName) (neighGroupName)
	stringstream comm(genComm);
	string boolean;
	comm>>boolean;
	int groupNN=0;
	Range nnRange;
	if (boolean=="equal"){
		comm>>groupNN;
	}
	else if (boolean=="inRange"){
		comm>>nnRange.x0>>nnRange.xf;
	}
	else if (boolean=="outOfRange"){
		comm>>nnRange.x0>>nnRange.xf;
	}
	else{
		child->NoGroupBoolean(boolean);
	}
	string sRName, neighGroupName;
	comm>>sRName>>neighGroupName;
	Range *sR=sRList.FindSpaceRange(sRName);
	Group *neighGroup=FindGroup(neighGroupName);
	parent->BuildNeighsInRange_fast(neighGroup, sR->x0, sR->xf);

	child->nAtoms=parent->nAtoms;
	child->AllocAtoms();
	int nCAtoms=0;
	for (int i=0; i<parent->nAtoms;i++){
		if (boolean=="equal"){
			if (parent->atom[i]->n_NN==groupNN){
				child->atom[nCAtoms]=parent->atom[i];
				nCAtoms++;
			}
		}
		else if (boolean=="inRange"||boolean=="outOfRange"){
			if (parent->atom[i]->n_NN>=nnRange.x0&&parent->atom[i]->n_NN<=nnRange.xf&&boolean=="inRange"){
				child->atom[nCAtoms]=parent->atom[i];
				nCAtoms++;
			}
			if ((parent->atom[i]->n_NN<nnRange.x0||parent->atom[i]->n_NN>nnRange.xf)&&boolean=="outOfRange"){
				child->atom[nCAtoms]=parent->atom[i];
				nCAtoms++;
			}
		}
		else{
			child->NoGroupBoolean(boolean);
		}
	}
	child->ShrinkGroup();
	child->p=parent->p;
	child->v=parent->v;
	cout<<"Created group with "<<child->nAtoms<<endl;

}
