/*
 * Timestep.cpp
 *
 *  Created on: Nov 18, 2010
 *      Author: ken
 */

#include "Timestep.h"

Timestep::Timestep() {
	step=-1;
	all.name="all";
	all.type="all";
	nUserGroups=0;
	userGroup=0;
	nextStep=0;
	nPosData=0;
	posData=0;
	nVelData=0;
	velData=0;
	firstPos=0;
	firstVel=0;
	curPos=0;
	curVel=0;
	nInPosActions=0;
	inPosAction=0;
}

Timestep::~Timestep() {
	if (nUserGroups!=0){
		for (int i=0;i<nUserGroups;i++){
			delete userGroup[i];
		}
		delete [] userGroup;
	}
	if (nPosData!=0){
		for (int i=0;i<nPosData;i++){
			delete posData[i];
		}
		delete [] posData;
	}
	if (nVelData!=0){
		for (int i=0;i<nVelData;i++){
			delete velData[i];
		}
		delete [] velData;
	}
	delete [] inPosAction;
}

void Timestep::InputPos(string posName)
{
	for (int i=0;i<nPosData;i++){
		if (posData[i]->name==posName||posName=="all"){
			if (posData[i]->p->fileType=="lammpstrj"){
				ReadPosInLammpstrjFile(posData[i], posData[i]->p);

				// Remap atoms using periodic boundary conditions
				for (int j=0; j<3; j++){
					if (posData[i]->p->boundary[j]=='p'){
						posData[i]->RemapPos(j,0,1);
					}
				}

				//Scale and translate position
				gsl_vector_view temp=gsl_vector_view_array(&posData[i]->p->sideL[0],3);
				posData[i]->ScalePos(&(temp.vector));
				gsl_vector *zero=gsl_vector_alloc(3);
				for (int j=0;j<3;j++){
					zero->data[j]=posData[i]->p->dim[j][0];
				}
				posData[i]->TransPos(zero);
				gsl_vector_free(zero);
			}
			else{
				cout<<"Position file type: "<<posData[i]->p->fileType<<" is not supported.\n";
				exit(0);
			}
			posData[i]->SortByIndex();
			DoInputPosActions(posData[i]);
		}
	}
	GenGroups();
}
void Timestep::DoInputPosActions(Group *inPos)
{
	for (int i=0;i<nInPosActions; i++){
		stringstream temp(inPosAction[i]);
		string name;
		temp>>name;
		if (inPos->name==name){
			string action;
			temp>>action;
			if (action=="scale"){
				double scaleFactor=1.0;
				temp >> scaleFactor;
				inPos->CalcCenter();
				gsl_vector *transVector=gsl_vector_alloc(3);
				gsl_vector_memcpy(transVector,inPos->p->center);
				gsl_vector_scale(transVector, -1.0);
				inPos->TransPos(transVector);
				gsl_vector *scaleVect=gsl_vector_alloc(3);
				gsl_vector_set_all(scaleVect, scaleFactor);
				inPos->ScalePos(scaleVect);
				inPos->TransPos(inPos->p->center);
				bool atomsIn=inPos->CheckAtomsInBox();
				if (!atomsIn){
					cout<<"Warning; new scaled position puts atoms outside box.\n --- If PBCs are used could result in overlapping atoms.\n";
				}
				gsl_vector_free(transVector);
				gsl_vector_free(scaleVect);
			}
			else{
				cout<< "Position data action: "<<action<<" not supported.\n";
				exit(0);
			}
		}
	}
}
void Timestep::InputVel(string velName)
{
	for (int i=0; i<nVelData; i++){
		if (velData[i]->name==velName||velName=="all"){
			if (velData[i]->v->fileType=="vel"){
				ReadVelInLammpsVelFile(velData[i], velData[i]->v);
			}
			else{
				cout<<"Velocity file type: "<<velData[i]->v->fileType<<" is not supported.\n";
				exit(0);
			}
			velData[i]->SortByIndex();
		}
	}
	GenGroups();
}
void Timestep::CompileAllVel()
{
	//Compile temporary list of all atoms from specified velocity files having that timestep
	if (nVelData){
		Group tempAll;
		for (int i=0;i<nVelData; i++){
			if (velData[i]->atom){
				tempAll.nAtoms+=velData[i]->nAtoms;
			}
		}
		if (tempAll.nAtoms){
			tempAll.AllocAtoms();
			int n=0;
			for (int i=0; i<nVelData; i++){
				if (velData[i]->atom){
					for (int j=0;j<velData[i]->nAtoms;j++){
						tempAll.atom[j+n]=velData[i]->atom[j];
					}
					n+=velData[i]->nAtoms;
				}
			}
			//Check if temporary velocity list matches any existing list of all atoms
			if (GroupsAreEqual(&tempAll, &all)||all.nAtoms==0){
				all.AllocVel();
				//Note: all contains a copy of the velocity data so 2x the data is currently saved in memory
				all.CopyAtomInfo(&tempAll);
				all.CopyAtomVel(&tempAll);
			}
			else{
				cout<<"Error in CompileAllPos, Groups not equal.\n";
				exit(0);
			}
		}
	}
}
void Timestep::CompileAllPos()
{
	//Compile temporary list of all atoms from any positions files which contain this timestep
	if (nPosData){
		Group tempAll;
		for (int i=0;i<nPosData; i++){
			if (posData[i]->atom){
				tempAll.nAtoms+=posData[i]->nAtoms;
			}
		}
		if (tempAll.nAtoms){
			tempAll.AllocAtoms();
			int n=0;
			for (int i=0; i<nPosData; i++){
				if (posData[i]->atom){
					for (int j=0;j<posData[i]->nAtoms;j++){
						tempAll.atom[j+n]=posData[i]->atom[j];
					}
					n+=posData[i]->nAtoms;
				}
			}
			//Check to make sure the atoms are the same atoms as those loaded from the any existing velocity info
			if (GroupsAreEqual(&tempAll, &all)||all.nAtoms==0){
				all.AllocPos();
				//Find Largest Box to fit all atoms
				for (int i=0;i<3;i++){
					if (nPosData){
						all.p->dim[i][0]=posData[0]->p->dim[i][0];
						all.p->dim[i][1]=posData[0]->p->dim[i][1];
						all.p->boundary[i]=posData[0]->p->boundary[i];
					}
					for (int j=1;j<nPosData;j++){
						all.p->dim[i][0]=getmin(all.p->dim[i][0], posData[j]->p->dim[i][0]);
						all.p->dim[i][1]=getmax(all.p->dim[i][1], posData[j]->p->dim[i][1]);
						if (all.p->boundary[i]!=posData[j]->p->boundary[i]){
							cout<<"Warning: Boundary conditions of simulation boxes are not equivalent.\n";
						}
					}
					all.p->sideL[i]=all.p->dim[i][1]-all.p->dim[i][0];
				}
				//Note: all is a copy of the position info so 2x the data are stored in memory
				all.CopyAtomInfo(&tempAll);
				all.CopyAtomPos(&tempAll);
			}
			else{
				cout<<"Error in CompileAllPos, Groups not equal.\n";
				exit(0);
			}
		}
	}
}
bool Timestep::GroupsAreEqual(Group *a, Group *b)
{
	if (a->nAtoms!=b->nAtoms){
		return false;
	}
	for (int i=0;i<a->nAtoms;i++){
		if (a->atom[i]->index!=b->atom[i]->index){
			cout<<"Atom indexes of input files are not equal.\n";
			return false;
		}
		if (a->atom[i]->type!=b->atom[i]->type){
			cout<<"Atom types of input files are not equal.\n";
			return false;
		}
	}
	return true;
}
void Timestep::GetMaxMinAtomCoords()
{
	for (int i=0;i<2;i++){
		for (int j=0;j<3;j++){
			all.p->dim[j][i]=gsl_vector_get(all.atom[0]->pos, j);
		}
	}
	for (int i=1;i<all.nAtoms;i++){
		for (int j=0;j<3;j++){
			double temp=gsl_vector_get(all.atom[i]->pos, j);
			all.p->dim[j][1]=getmax(temp, all.p->dim[j][1]);
			all.p->dim[j][0]=getmin(temp, all.p->dim[j][0]);
		}
	}
}

void Timestep::BuildNearestNeighs(double Rc)
{
	double RcSq=Rc*Rc;
	double searchBound[3][2];
	int tempn_NN=0;
	Atom **tempNN=new Atom* [all.nAtoms];
	//TODO consider removing zeroing and just copying nPointers
	//Create temporary NN list which is then copied over to each atom
	for (int i=0;i<all.nAtoms;i++){
		tempNN[i]=0;
	}
	for (int i=0;i<all.nAtoms;i++){
		for (int j=0; j<tempn_NN;j++){
			tempNN[j]=0;
		}
		tempn_NN=0;
		for (int dim=0;dim<3;dim++){
			double x=gsl_vector_get(all.atom[i]->pos, dim);
			searchBound[dim][0]=x-Rc;
			searchBound[dim][1]=x+Rc;
		}
		//Check NN lists of already calculated atoms
		for (int j=0;j<i;j++){
			for (int k=0;k<all.atom[j]->n_NN;k++){
				if (all.atom[j]->NN[k]==all.atom[i]){
					tempNN[tempn_NN]=all.atom[j];
					tempn_NN++;
				}
			}
		}
		//Find new NN
		for (int j=i+1; j<all.nAtoms;j++){
			bool inBounds=true;
			// To diminish unnecessary distance calcs, create box around each atom and use inequalities.
			for (int dim=0;dim<3;dim++){
				double x=gsl_vector_get(all.atom[j]->pos,dim);
				if (x<searchBound[dim][0]||x>searchBound[dim][1]){
					inBounds=false;
				}
			}
			if (inBounds){
				//TODO Consider PBC in distance calc
				gsl_vector *temp=gsl_vector_alloc(3);
				double dSq=0;
				gsl_vector_memcpy(temp, all.atom[i]->pos);
				gsl_vector_sub(temp, all.atom[j]->pos);
				gsl_blas_ddot(temp, temp, &dSq);
				if (dSq<=RcSq){
					tempNN[tempn_NN]=all.atom[j];
					tempn_NN++;
				}
				gsl_vector_free(temp);
			}
		}
		//Copy NN list to atoms
		all.atom[i]->n_NN=tempn_NN;
		AllocAtomNN(all.atom[i]);
		for(int j=0;j<tempn_NN;j++){
			all.atom[i]->NN[j]=tempNN[j];
		}
	}
	delete [] tempNN;
}
void Timestep::BuildNeighsInRange(double R0, double Rc)
{
	double RcSq=Rc*Rc, R0Sq=R0*R0;
	double innerBound[3][2], outerBound[3][2];
	int tempn_NN=0;
	Atom **tempNN=new Atom* [all.nAtoms];
	//TODO consider removing zeroing and just copying nPointers
	//Create temporary NN list which is then copied over to each atom
	for (int i=0;i<all.nAtoms;i++){
		tempNN[i]=0;
	}
	for (int i=0;i<all.nAtoms;i++){
		for (int j=0; j<tempn_NN;j++){
			tempNN[j]=0;
		}
		tempn_NN=0;
		// outerBound is cube inscribing a sphere of radius Rc
		// innerBound is cube inscribed in sphere of radius R0
		for (int dim=0;dim<3;dim++){
			double x=gsl_vector_get(all.atom[i]->pos, dim);
			outerBound[dim][0]=x-Rc;
			outerBound[dim][1]=x+Rc;
			innerBound[dim][0]=x-(R0/sqrt(3.0));
			innerBound[dim][0]=x+(R0/sqrt(3.0));
		}
		//Check NN lists of already calculated atoms
		for (int j=0;j<i;j++){
			for (int k=0;k<all.atom[j]->n_NN;k++){
				if (all.atom[j]->NN[k]==all.atom[i]){
					tempNN[tempn_NN]=all.atom[j];
					tempn_NN++;
				}
			}
		}
		//Find new NN
		for (int j=i+1; j<all.nAtoms;j++){
			bool inBounds=true;
			// To diminish unnecessary distance calcs, create box around each atom and use inequalities.
			for (int dim=0;dim<3;dim++){
				double x=gsl_vector_get(all.atom[j]->pos,dim);
				if (x<outerBound[dim][0]||x>outerBound[dim][1]||(x>innerBound[dim][0]&&x<innerBound[dim][1])){
					inBounds=false;
				}
			}
			if (inBounds){
				//TODO Consider PBC in distance calc
				gsl_vector *temp=gsl_vector_alloc(3);
				double dSq=0;
				gsl_vector_memcpy(temp, all.atom[i]->pos);
				gsl_vector_sub(temp, all.atom[j]->pos);
				gsl_blas_ddot(temp, temp, &dSq);
				if (dSq<=RcSq&&dSq>=R0Sq){
					tempNN[tempn_NN]=all.atom[j];
					tempn_NN++;
				}
				gsl_vector_free(temp);
			}
		}
		//Copy NN list to atoms
		all.atom[i]->n_NN=tempn_NN;
		AllocAtomNN(all.atom[i]);
		for(int j=0;j<tempn_NN;j++){
			all.atom[i]->NN[j]=tempNN[j];
		}
	}
	delete [] tempNN;
}
Group * Timestep::FindGroup(string name)
{
	if (name=="all"){
		return &all;
	}
	for (int i=0;i<nUserGroups;i++){
		if (userGroup[i]->name==name){
			return userGroup[i];
		}
	}
	cout<<"Cannot find group " <<name<<" in timestep "<<step<<endl;
	exit (0);
}
Group* Timestep::FindPos(string name)
{
	for (int i=0; i<nPosData;i++){
		if (posData[i]->name==name){
			return posData[i];
		}
	}
	return 0;
}
//TODO Consider sorting pos and vel array to make sure that they correspond
void Timestep::CreatePosArray()
{
	try{
		if (posData){
			delete [] posData;
		}
		posData=new Group *[nPosData];
	}
	catch(exception &e){
		exit(1);
	}
	Group *curGroup=firstPos;
	for (int i=0;i<nPosData;i++){
		posData[i]=curGroup;
		curGroup=curGroup->next;
	}
}
void Timestep::CreateVelArray()
{
	try{
		if (velData){
			delete [] velData;
		}
		velData=new Group *[nVelData];
	}
	catch(exception &e){
		exit(1);
	}
	Group *curGroup=firstVel;
	for (int i=0;i<nVelData;i++){
		velData[i]=curGroup;
		curGroup=curGroup->next;
	}
}
void Timestep::ReadVelInLammpsVelFile(Group *tempAtoms, Velocity *tempV)
{
	string fileName=tempV->file;
	ifstream infile;
	infile.open(fileName.c_str());
	if (infile.is_open()){
		infile.seekg(tempV->dataPos);
		string tempstring, item;
		bool stop=false;
		while(!infile.eof()&&!stop){
			infile >> tempstring >>item;
			if(item=="NUMBER"){
				infile.ignore(100,'\n');
				infile >> tempAtoms->nAtoms;
				//allocate mem for positions and other info
				tempAtoms->AllocAtoms();
				tempAtoms->AllocAtomVel();
			}
			else if (item=="BOX"){
				infile.ignore(100,'\n');
				for (int i=0;i<3;i++){
					infile.ignore(100,'\n');
				}
			}
			//Input atom positions
			else if (item=="ATOMS"){
				infile.ignore(100,'\n');
				for (int i=0; i<tempAtoms->nAtoms; i++){
					infile >>tempAtoms->atom[i]->index>>tempAtoms->atom[i]->type>>tempAtoms->atom[i]->vel->data[0]>> tempAtoms->atom[i]->vel->data[1] >>tempAtoms->atom[i]->vel->data[2];
				}
				stop=true;
			}
			else{
				cout<<"Unsupported keyword in lammps velocity file: "<<item<<endl;
				exit(0);
			}
		}
		infile.close();
	}
	else {
		printf("Could not open velocity file: %s \n", fileName.c_str());
		exit(0);
	}
}
void Timestep::ReadPosInLammpstrjFile(Group *tempAtoms, Position *tempP)
{
	if (!tempAtoms){
		cout<<"Error in pos file read, group of atoms uninitialized."<<endl;
		exit(1);
	}
	ifstream infile;
	infile.open(tempP->file.c_str());
	if (infile.is_open()){
		infile.seekg(tempP->dataPos);
		string tempstring, item;
		bool stop=false;
		while (!infile.eof()&&!stop){
			infile >>tempstring >>item;
			if(item=="NUMBER"){
				infile.ignore(100,'\n');
				infile >> tempAtoms->nAtoms;
				//allocate mem for positions and other info
				tempAtoms->AllocAtoms();
				tempAtoms->AllocAtomPos();
			}
			else if (item=="BOX"){
				infile.ignore(100,'\n');
				for (int i=0;i<3;i++){
					for(int j=0;j<2;j++){
						infile>>tempP->dim[i][j];
					}
					tempP->sideL[i]=tempP->dim[i][1]-tempP->dim[i][0];
				}
			}
			//Input atom positions
			else if (item=="ATOMS"){
				infile.ignore(100,'\n');
				for (int i=0; i<tempAtoms->nAtoms; i++){
					infile >>tempAtoms->atom[i]->index>>tempAtoms->atom[i]->type>>tempAtoms->atom[i]->pos->data[0]>> tempAtoms->atom[i]->pos->data[1] >>tempAtoms->atom[i]->pos->data[2];
				}
				stop=true;
			}
			else{
				cout<<"Unsupported keyword in lammpstrj file: "<<item<<endl;
				exit(0);
			}
		}
		infile.close();
	}
	else {
		printf("Could not open file: %s \n", tempP->file.c_str());
		exit(0);
	}
}


void Timestep::Clean()
{
	/*
	if (subDomain){
		for(int i=0;i<nSubDomX;i++){
			for (int j=0; j<nSubDomY; j++){
				for (int k=0;k<nSubDomZ;k++){
					delete [] subDomain[i][j][k].atom;
				}
				nSubDomZ=0;
				delete [] subDomain[i][j];
			}
			nSubDomY=0;
			delete [] subDomain[i];
		}
		nSubDomX=0;
		delete [] subDomain;
	}
	*/
	//Deletes atom information in all but does not destroy the group all

	for (int i=0;i<all.nAtoms;i++){
		if (all.atom[i]->pos){
			gsl_vector_free(all.atom[i]->pos);
		}
		if (all.atom[i]->vel){
			gsl_vector_free(all.atom[i]->vel);
		}
		if (all.atom[i]->n_NN){
			delete [] all.atom[i]->NN;
		}
		delete all.atom[i];
	}
	delete [] all.atom;
	all.atom=0;
	all.nAtoms=0;
	// Deletes atom information from any loaded position file
	for (int i=0;i<nPosData;i++){
		if (posData[i]->atom){
			for(int j=0;j<posData[i]->nAtoms;j++){
				if (posData[i]->atom[j]->pos){
					gsl_vector_free(posData[i]->atom[j]->pos);
				}
				delete posData[i]->atom[j];
			}
			delete [] posData[i]->atom;
			posData[i]->atom=0;
			posData[i]->nAtoms=0;
		}
	}
	//Deletes atom information from any loaded velocity file
	for (int i=0;i<nVelData;i++){
		if (velData[i]->atom){
			for(int j=0;j<velData[i]->nAtoms;j++){
				if (velData[i]->atom[j]->vel){
					gsl_vector_free(velData[i]->atom[j]->vel);
				}
				delete velData[i]->atom[j];
			}
			delete [] velData[i]->atom;
			velData[i]->atom=0;
			velData[i]->nAtoms=0;
		}
	}
}
// TODO Function will not work and needs work!
void Timestep::CalcLatticeParameter(Lattice *latt)
{
	Group centerAtoms;
	centerAtoms.nAtoms=0;
	all.CalcCenter();
	all.CalcPosMags(all.p->center);
	double r_c=0;
	for (int i=0;i<3;i++){
		double x=0;
		for (int j=0;j<3;j++){
			x+=all.p->lattice.vectors[j][i];
		}
		r_c+=x*x;
	}
	r_c=sqrt(r_c)/2.0;
	for(int i=0;i<all.nAtoms;i++){
		if(all.atom[i]->r<r_c){
			centerAtoms.nAtoms++;
		}
	}
	centerAtoms.AllocAtoms();
	int count=0;
	for (int i=0;i<all.nAtoms;i++){
		if(all.atom[i]->r<r_c){
			centerAtoms.atom[count]=all.atom[i];
			count++;
		}
	}
	//calculating full RDF and then averaging values near
	//expected nearest neighbors to determine new unit cell parameter.

	all.p->lattice.CopyLattice(latt);
	double nn, nn2, a_nn, a_nn2, del;
	if (latt->name=="faceCentered"&&latt->symmetry=="cubic"){
		nn=latt->a/sqrt(2.0);
		nn2=latt->a;
		del=nn*.15;
	}
	else{
		cout<<"Lattice: "<<latt->name<<" and symmetry: "<<latt->symmetry<<" not supported in calc of new lattice parameter.\n";
		exit(0);
	}
	//TODO Will not work because bounds not set
	centerAtoms.CalcRDF(DPREC/10.0);
	double avgNN=0, avgNN2=0;
	double totNN=0, totNN2=0;
	//Calc nn by averaging
	for (int i=(int)((nn-del)/centerAtoms.p->rdf->step);i<(int)((nn+del)/centerAtoms.p->rdf->step);i++){
		avgNN+=centerAtoms.p->rdf->x[i]*centerAtoms.p->rdf->y[i];
		totNN+=centerAtoms.p->rdf->y[i];
	}
	for (int i=(int)((nn2-del)/centerAtoms.p->rdf->step);i<(int)((nn2+del)/centerAtoms.p->rdf->step);i++){
		avgNN2+=centerAtoms.p->rdf->x[i]*centerAtoms.p->rdf->y[i];
		totNN2+=centerAtoms.p->rdf->y[i];
	}
	if (latt->name=="faceCentered"&&latt->symmetry=="cubic"){
		a_nn=avgNN*sqrt(2.0)/(double) totNN;
		a_nn2=avgNN2/(double)totNN2;
		latt->a=latt->b=latt->c=(a_nn+a_nn2)/2.0;
	}
	latt->InitLatticeVect();
	CleanDistrib(centerAtoms.p->rdf);
}
void Timestep::DeletePosData(string name)
{
	for (int i=0;i<nPosData;i++){
		if (posData[i]->name==name){
			if (i==0){
				firstPos=posData[i]->next;
			}
			else if (i==nPosData-1){
				posData[i-1]->next=0;
			}
			else{
				posData[i-1]->next=posData[i+1];
			}
			delete posData[i];
			nPosData--;
		}
	}
	CreatePosArray();
}
void Timestep::CreatePosInputAction(string command)
{
	//Add a command to the list to be executed when inputing a position
	//Create a copy of existing commands and swap pointers at end.
	string *newCommands=new string [nInPosActions+1];
	for (int i=0;i<nInPosActions;i++){
		newCommands[i]=inPosAction[i];
	}
	newCommands[nInPosActions]=command;
	nInPosActions++;
	delete [] inPosAction;
	inPosAction=newCommands;
}
