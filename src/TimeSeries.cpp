/*
 * TimeSeries.cpp
 *
 *  Created on: Jan 11, 2010
 *      Author: Ken
 */

#include "TimeSeries.h"
TimeSeries::TimeSeries() {
	nTimeSteps=0;
	name="";
	lattice=0;

	relaxedLatt=0;
	staticDisp=0;
	dynDisp=0;
	userDisp=0;
	nDynDisp=0;
	avgDynDisp=0;
	nStrTen=0;
	F=0;

	firstStep=0;
	step=0;
	velCorrFn=0;
	avg=0;
	msd_t=0;

	nTimeRange=0;
	tRange=0;
}

TimeSeries::~TimeSeries()
{
	if (F){
		for (int i=0;i<nStrTen;i++){
			CleanAtomTensor(&(F[i]));
		}
		delete [] F;
	}
	if (staticDisp){
		delete staticDisp;
	}
	if (userDisp){
		delete userDisp;
	}
	if (dynDisp){
		for(int i=0; i<nDynDisp; i++){
			delete dynDisp[i];
		}
		delete [] dynDisp;
	}
	if (msd_t){
		CleanDistrib(msd_t);
		delete msd_t;
	}

	if (lattice){
		delete lattice;
	}
	if(relaxedLatt){
		delete relaxedLatt;
	}
	if (step){
		for (int i=0;i<nTimeSteps;i++){
			delete step[i];
		}
		delete [] step;
	}
	if (velCorrFn){
		CleanDistrib(velCorrFn);
		delete velCorrFn;
	}
	if (avg){
		delete avg;
	}
	if (tRange){
		for (int i=0;i<nTimeRange;i++){
			delete tRange[i];
		}
		if (nTimeRange>1){
			delete [] tRange;
		}
		else{
			delete tRange;
		}

	}
	// Clean up the static variables of list of space ranges
	sRList.Clean();
}

void TimeSeries::FindPosInLammpstrjFile(string posGroupName, string file, char *boundary)
{
	string fileName=file;
	Timestep *curTime=firstStep, *prevTime=0;
	string tempstring, item;
	ifstream infile;
	infile.open(fileName.c_str());
	if (infile.is_open()){
		printf("Lammpstrj file %s opened\n", fileName.c_str());
		while (!infile.eof())
		{
			infile >> tempstring;
			if(!infile.eof()){
				if (tempstring=="ITEM:"){
					infile>>item;
					//Find locations of timesteps in file
					if (item=="TIMESTEP"){
						int tempTime;
						infile >>tempTime;
						InsertStepInSeq(&prevTime, &curTime, tempTime);
						curTime->step=tempTime;
						//add position to list of datafiles for timestep
						Group *newPos=new Group;
						if (curTime->curPos){
							curTime->curPos->next=newPos;
						}
						curTime->curPos=newPos;
						curTime->nPosData++;
						if (!curTime->firstPos){
							curTime->firstPos=curTime->curPos;
						}
						//
						curTime->curPos->AllocPos();
						curTime->curPos->name=posGroupName;
						curTime->curPos->type="pos";
						curTime->curPos->p->file=fileName;
						curTime->curPos->p->dataPos=infile.tellg();
						curTime->curPos->p->fileType="lammpstrj";
						for (int i=0;i<3;i++){
							curTime->curPos->p->boundary[i]=boundary[i];
						}
						curTime->CreatePosArray();
						cout<<"Found position timeStep: "<<curTime->step<<" in file: "<<file<<endl;
						prevTime=curTime;
						curTime=curTime->nextStep;
					}
					else{
						infile.ignore(100,'\n');
					}
				}
				else{
					infile.ignore(100,'\n');
				}
			}
		}
		infile.close();
	}
	else {
		printf("Could not open file: %s \n", fileName.c_str());
		exit(0);
	}
}
void TimeSeries::FindVelInLammpsVelFile(string velGroupName, string file)
{
	string fileName=file;
	Timestep *curTime=firstStep, *prevTime=0;
	string tempstring, item;
	ifstream infile;
	infile.open(fileName.c_str());
	if (infile.is_open()){
		printf("Lammps .vel file %s opened\n", fileName.c_str());
		while (!infile.eof())
		{
			infile >> tempstring;
			if(!infile.eof()){
				if (tempstring=="ITEM:"){
					infile>>item;
					//Find locations of timesteps in file
					if (item=="TIMESTEP"){
						int tempTime;
						infile >>tempTime;
						InsertStepInSeq(&prevTime, &curTime, tempTime);
						curTime->step=tempTime;
						//add velocity to list of velocity files for timestep
						Group *newGroup=new Group;
						if (curTime->curVel){
							curTime->curVel->next=newGroup;
						}
						curTime->curVel=newGroup;
						curTime->nVelData++;
						if(!curTime->firstVel){
							curTime->firstVel=curTime->curVel;
						}
						//
						curTime->curVel->AllocVel();
						curTime->curVel->name=velGroupName;
						curTime->curVel->type="vel";
						curTime->curVel->v->file=fileName;
						curTime->curVel->v->fileType="vel";
						curTime->curVel->v->dataPos=infile.tellg();
						curTime->CreateVelArray();
						cout<<"Found velocity timeStep: "<<curTime->step<<" in file: "<<file<<endl;
						prevTime=curTime;
						curTime=curTime->nextStep;
					}
					else{
						infile.ignore(100,'\n');
					}
				}
				else{
					infile.ignore(100,'\n');
				}
			}
		}
		infile.close();
	}
	else {
		printf("Could not open file: %s \n", fileName.c_str());
		exit(0);
	}
}
void TimeSeries::InsertStepInSeq(Timestep **prev, Timestep **cur, int t)
{
	if (*cur==0){
		*cur=new Timestep();
		nTimeSteps++;
		if (*prev!=0){
			(*prev)->nextStep=*cur;
		}
		if (firstStep==0){
			firstStep=*cur;
		}
	}
	else{
		//If t>curTime move to next step and recall function
		if(t>(*cur)->step){
			*prev=*cur;
			*cur=(*cur)->nextStep;
			InsertStepInSeq(prev, cur, t);
		}
		//If t<curTime insert new timestep
		else if(t<(*cur)->step){
			Timestep *nextStep=*cur;
			*cur=new Timestep;
			nTimeSteps++;
			(*cur)->nextStep=nextStep;
			if (*prev!=0){
				(*prev)->nextStep=*cur;
			}
			else{
				firstStep=*cur;
			}
		}
		//If t=curTime return pointers
		else{
			return;
		}
	}
}
void TimeSeries::CreateSeriesArray()
{
	try{
		if (step){
			delete [] step;
		}
		step=new Timestep *[nTimeSteps];
	}
	catch(exception &e){
		cout<<e.what()<<endl;
		exit(1);
	}
	Timestep *cur=firstStep;
	for (int i=0;i<nTimeSteps;i++){
		step[i]=cur;
		cur=cur->nextStep;
	}
}
int TimeSeries::GetSmallConstTimeStep(int start, int finish, string type)
{
	int del=0;
	int t1=-1;
	int t2=-1;
	for(int i=0;i<nTimeSteps;i++){
		t2=t1;
		if (type=="vel"){
			if (step[i]->firstVel){
				t1=step[i]->step;
			}
		}
		else if (type=="pos"){
			if (step[i]->firstPos){
				t1=step[i]->step;
			}
		}
		else if (type=="all"){
			t1=step[i]->step;
		}
		else{
			cout<<"Type not supported in GetSmallConstTimestep: "<<type<<"\n";
			exit(0);
		}
		if (t2>=start && t1<=finish){
			if (del!=0){
				if (del!=(t1-t2)){
					cout<<"Error: no constant step in the range: "<<start<<" - "<<finish<<endl;
					exit(0);
				}
			}
			else{
				del=t1-t2;
			}
		}
	}
	cout<<"Found min timeStep of :" <<del<<endl;
	return del;
}

/*
void TimeSeries::MakeDistribs()
{
	if(calcDynDisp){
		if (calcAvgDynDisp){
			AllocAvgDynDispDistrib();
			CalcAvgDynDispDistrib();
		}
	}
}
*/
void TimeSeries::CalcDisplace(Displace **d, string _specifyMode, string _groupName, bool scaled, int t0, int tf)
{
	try{
		if (!*d){
			*d=new Displace;
		}
	}
	catch(exception &e){
		cout<<e.what()<<endl;
	}
	(*d)->initial=t0;
	(*d)->final=tf;
	(*d)->specifyMode=_specifyMode;
	(*d)->groupName=_groupName;
	(*d)->scaled=scaled;
	AssignDispPositions(*d);
	(*d)->Allocate();
	//Calc displacement from old lattice parameter.
	if ((*d)->scaled){
		//TODO write support for calculated lattice parameter
		if (lattice->name=="faceCentered"&&lattice->symmetry=="cubic"){
			(*d)->CalcDisplacements(1.0/lattice->a, 1.0/relaxedLatt->a);
		}
	}
	else {
		(*d)->CalcDisplacements(1.0, 1.0);
	}
}
/*
void TimeSeries::CopyPosition(Position *from, Position *to)
{
	to->nAtoms=from->nAtoms;
	to->Allocate();
	for(int i=0;i<3;i++){
		for(int j=0;j<2;j++){
			to->dim[i][j]=from->dim[i][j];
		}
		to->center[i]=from->center[i];
	}
	for(int i=0;i<to->nAtoms;i++){
		to->CopyAtom(&from->atom[i],&to->atom[i]);
	}
}
*/

//TODO make sure avg is always allocated before using this function!
void TimeSeries::AssignDispPositions(Displace *d)
{
	if (d->initial==-1){
		d->pos_i=avg;
	}
	if (d->final==-1){
		d->pos_f=avg;
	}
	if (d->specifyMode=="timeStep"){
		for (int j=0;j<nTimeSteps;j++){
			if(d->initial==step[j]->step){
				d->pos_i=step[j]->FindGroup(d->groupName);
			}
			if(d->final==step[j]->step){
				d->pos_f=step[j]->FindGroup(d->groupName);
			}
		}
		if(d->pos_i==0||d->pos_f==0){
			cout<<"Could not find timeSteps: "<<d->initial<< " and "<<d->final<<endl;
			exit(0);
		}
	}
	else if (d->specifyMode=="frame"){
		if (d->initial>0&&d->initial<nTimeSteps+1){
			d->pos_i=step[d->initial-1]->FindGroup(d->groupName);
		}
		else{
			cout<<"Could not find frame: " <<d->initial<<endl;
			exit(0);
		}
		if (d->final>0&&d->final<nTimeSteps+1){
			d->pos_f=step[d->final-1]->FindGroup(d->groupName);
		}
		else{
			cout<<"Could not find frame: " <<d->final<<endl;
			exit(0);
		}
	}
	else{
		cout<<"Displacement specifyMode unsupported: "<<d->specifyMode<<endl;
	}
}
//TODO Only works with input timestep, not frames!
void TimeSeries::CalcStaticDisp(string groupName, bool scaled, int ref, int t0, int tf)
{
	cout<<"\nCalculating Static Displacement... ";
	CalcAvgPos(groupName, t0, tf);
	try{
		staticDisp=new Displace;
	}
	catch(exception &e){
		cout<<"Exception in TimeSeries::CalcStaticDisp: "<<e.what()<<endl;
		exit(1);
	}
	Timestep *t=FindTimestep(ref);
	if (t==0){
		exit(0);
	}
	if(t->firstPos!=0){
		t->InputPos("all");
		CalcDisplace(&staticDisp, "timeStep", groupName, scaled, ref, -1);
		t->Clean();
	}
	cout<<"done.\n";
}
void TimeSeries::CalcDynamicDisp(string groupName, bool scaled, int t0, int tf, int dt)
{
	cout<<"\nCalculating Dynamic Displacement ...";
	CalcAvgPos(groupName, t0, tf);
	nDynDisp=GetNumPosInRange(t0,tf,dt);
	try{
		dynDisp=new Displace* [nDynDisp];
		for (int i=0;i<nDynDisp;i++){
			dynDisp[i]=0;
		}
	}
	catch(exception &e){
		cout<<"Exception in TimeSeries::CalcDynamicDisp: "<<e.what()<<endl;
		exit(1);
	}
	int count=0;
	for(int i=0;i<nTimeSteps;i++){
		if (step[i]->step>=t0&&step[i]->step<=tf&&step[i]->firstPos!=0){
			step[i]->InputPos("all");
			CalcDisplace(&(dynDisp[count]), "timeStep",groupName, scaled, -1, step[i]->step);
			count++;
			step[i]->Clean();
			t0+=dt;
		}
	}
	cout<<"done.\n";
}
void TimeSeries::CalcMSD_t(string group,int ref, int t0, int tf )
{
	int n=GetNumPosInRange(t0,tf);
	if (n>0){
		cout<<"\nCalculating MSD(t)... ";
		AllocDistrib(&msd_t, "MSD(t)", "timestep", "MSD", 0, n, 1.0);
		Displace *tempDisp=0;
		Timestep *t=0;
		if (ref==-1){
			CalcAvgPos(group, t0, tf);
		}
		else{
			t=FindTimestep(ref);
			if (t==0){
				exit(0);
			}
			t->InputPos("all");
		}
		int count=0;
		for (int i=0;i<nTimeSteps;i++){
			if (step[i]->step>=t0&& step[i]->step<=tf&&step[i]->firstPos!=0){
				step[i]->InputPos("all");
				CalcDisplace(&tempDisp, "timeStep", group, false, ref, step[i]->step);
				tempDisp->CalcAvgMSD();
				msd_t->x[count]=step[i]->step;
				msd_t->y[count]=tempDisp->msd;
				count++;
				step[i]->Clean();
			}
		}
		if (t){
			t->Clean();
		}
		delete tempDisp;
		cout<<"done\n";
	}
}

/*
void TimeSeries::AllocAvgDynDispDistrib()
{
	try{
		if (avgDynDisp) delete avgDynDisp;
		avgDynDisp=new Displace;
	}
	catch(exception &e){
		cout<<"Exception in TimeSeries::CalcAvgDynDisp: "<<e.what()<<endl;
		exit(1);
	}
	avgDynDisp->AllocateMSD(GetMaxValOfXAxis(dynDisp[0].pMSD->name)+1);
	avgDynDisp->AllocateRadMSD(GetMaxValOfXAxis(dynDisp[0].radMSD->name)+1);
	avgDynDisp->AllocateRadDisp(GetMaxValOfXAxis(dynDisp[0].radDisp->name)+1);
	avgDynDisp->scaled=dynDisp[0].scaled;
	avgDynDisp->initial=-2;
	avgDynDisp->final=-2;

}
double TimeSeries::GetMaxValOfXAxis(string distribName)
{
	double max = 0;
	for (int i=0;i<nRelaxedPos;i++){
		if (distribName==dynDisp[i].pMSD->name){
			max=getmax(max,dynDisp[i].pMSD->maxX);
		}
		else if (distribName==dynDisp[i].radMSD->name){
			max=getmax(max, dynDisp[i].radMSD->maxX);
		}
		else if (distribName==dynDisp[i].radDisp->name){
			max=getmax(max, dynDisp[i].radDisp->maxX);
		}
		else {
			cout<<"Error in GetMaxValOfXAxis: Cannot find distribution named: "<<distribName<<endl;
			exit(0);
		}
	}
	return max;
}
void TimeSeries::CalcAvgDynDispDistrib()
{
	for(int i=0;i<nRelaxedPos;i++){
		for (int j=0;j<dynDisp[i].pMSD->n;j++){
			avgDynDisp->pMSD->y[(int)(dynDisp[i].pMSD->x[j]/avgDynDisp->pMSD->step)]+=dynDisp[i].pMSD->y[j];
		}
		for (int j=0;j<dynDisp[i].radMSD->n;j++){
			avgDynDisp->radMSD->y[(int)(dynDisp[i].radMSD->x[j]/avgDynDisp->radMSD->step)]+=dynDisp[i].radMSD->y[j];
		}
		for (int j=0;j<dynDisp[i].radDisp->n;j++){
			avgDynDisp->radDisp->y[(int)(dynDisp[i].radDisp->x[j]/avgDynDisp->radDisp->step)]+=dynDisp[i].radDisp->y[j];
		}
	}
	for (int i=0;i<avgDynDisp->pMSD->n;i++){
		avgDynDisp->pMSD->y[i]/=(double)nRelaxedPos;
	}
	for (int i=0;i<avgDynDisp->radMSD->n;i++){
		avgDynDisp->radMSD->y[i]/=(double)nRelaxedPos;
	}
	for (int i=0;i<avgDynDisp->radDisp->n;i++){
		avgDynDisp->radDisp->y[i]/=(double)nRelaxedPos;
	}
}
*/
void TimeSeries::CalcInstantLattice(int t0)
{
	double nn=lattice->GetNNDist(), nn2=lattice->GetNN2Dist();
	for (int i=0;i<nTimeSteps;i++){
		if (step[i]->step>=t0){
			lattice->CopyLattice(&(step[i]->all.p->lattice));
			step[i]->BuildNearestNeighs((nn+nn2)/2.0);
		}
	}
}

void TimeSeries::CalcAtomicStrainTensors(string group, int t0, int tf, int dt, double Rc, string weightType)
{
	Timestep *prevPos=0, *curPos=0;
	int prevTime=0, curTime=t0;
	cout<<"\nCalculating Atomic Strain ...";
	if (F){
		for (int i=0;i<nStrTen;i++){
			CleanAtomTensor(&(F[i]));
		}
		delete [] F;
	}
	nStrTen=GetNumPosInRange(t0, tf, dt)-1;
	try{
		F=new AtomStrain [nStrTen];
	}
	catch(exception &e){

	}
	int count=0;
	for (int i=0;i<nTimeSteps; i++){
		if (step[i]->step>=curTime&&step[i]->step<=tf&&step[i]->firstPos!=0){
			curPos=step[i];
			curPos->InputPos("all");

			if (prevPos){
				Group *g_i=prevPos->FindGroup(group), *g_f=curPos->FindGroup(group);
				if (g_i->nAtoms==g_f->nAtoms){
					F[count].nA=g_i->nAtoms;
					F[count].index=new int [F[count].nA];
					for (int j=0;j<g_i->nAtoms;j++){
						F[count].index[j]=g_i->atom[j]->index;
					}
				}
				else{
					cout<< "Cannot calculate atomic strain between t0: "<<t0<<" and tf: "<<tf<<" because not equal number of atoms.\n";
					exit(0);
				}
				F[count].t0=prevPos->step;
				F[count].tf=curPos->step;
				prevPos->BuildNearestNeighs(Rc);
				CalcStrainTensor(g_i, g_f, Rc, weightType, &F[count]);
				prevPos->Clean();
				count++;
			}
			prevPos=curPos;
			prevTime=step[i]->step;
			curTime+=dt;
		}
	}
	if(curPos){
		curPos->Clean();
	}
	cout<<"done.\n";
}

/* Pseudo-code for strain tensor calc
 * 1 find NN of atoms at t_i
 * 2 Calc NN vectors from t_i and same index for t_f
 * *** Ideally need to calc strain by comparing "lattice position", not atom pos with same index (future problem)
 * 3 Calc weights for NN
 * 4 Calc Local Deformation Gradiant Tensor
 * 5 Calculate Average Deformation Gradient Tensor
 * 6 Decompose Avg Def Grad Tensor into Rotation (R) and Strain(S)
 * 7 Apply R^-1 to particle before calc of displacements
 */
void TimeSeries::CalcStrainTensor(Group *g_i, Group *g_f, double Rc, string weightType, AtomStrain *atomF)
{
	try{
		atomF->s=new gsl_matrix *[atomF->nA];
	}
	catch (exception &e){
	}
	atomF->g_i=g_i;

	// Initialize NN Vectors
	//TODO Go back to refine here
	//lattice->CopyLattice(&(g_i->p->lattice));
	for (int i=0;i<g_i->nAtoms;i++){
		//Find atoms of same index as Group i in Group f
		Atom *g_f_atom_i=g_f->FindAtom(g_i->atom[i]->index);
		gsl_vector **NNti, **NNtf, *weight;
		try{
			NNti=new gsl_vector* [g_i->atom[i]->n_NN];
			NNtf=new gsl_vector* [g_i->atom[i]->n_NN];
			for (int j=0;j<g_i->atom[i]->n_NN;j++){
				NNti[j]=gsl_vector_calloc(3);
				NNtf[j]=gsl_vector_calloc(3);
			}
			weight=gsl_vector_calloc(g_i->atom[i]->n_NN);
			atomF->s[i]=gsl_matrix_calloc(3,3);
		}
		catch(exception &e){
			cout<<"Exception in CalcStrainTensor: "<<e.what()<<endl;
			exit(1);
		}
		for (int j=0;j<g_i->atom[i]->n_NN;j++){
			Atom *g_f_NNj=g_f->FindAtom(g_i->atom[i]->NN[j]->index);
			gsl_vector_memcpy(NNti[j], g_i->atom[i]->NN[j]->pos);
			gsl_vector_sub(NNti[j],g_i->atom[i]->pos);
			gsl_vector_memcpy(NNtf[j], g_f_NNj->pos);
			gsl_vector_sub(NNtf[j],g_f_atom_i->pos);
		}
		//Calculate weights
		if (weightType=="HeavysideStep"){
			ApplyHeavysideStepWeights(weight, NNti, g_i->atom[i]->n_NN, Rc);
		}
		else{
			cout<<"Weight type: "<<weightType<<" is not supported in CalcStrainTensor.\n";
			exit(0);
		}
		// Calc Local Strain Gradent
		CalcLocalDefGradTensor(NNti, NNtf, weight, g_i->atom[i]->n_NN, atomF->s[i]);
		//Cleanup
		gsl_vector_free(weight);
		for (int j=0;j<g_i->atom[i]->n_NN;j++){
			gsl_vector_free(NNti[j]);
			gsl_vector_free(NNtf[j]);
		}
		delete [] NNti;
		delete [] NNtf;
	}
}

//Method of calculating the Local Deformaiton Gradient Tensor from the nearest neighbor environment of two different
//	configurations which is described in Gullet, Modelling Sim. Mat. Sci Eng. 16 (2008) 015001
//As written the routine takes in an array of NN distance vectors for each configuration and a vector of weights
//  *** It is assumed that the arrays are already ordered and that each NN is identified.
// Output is the "optimal" local deformation tensor, F, found by least squares of distance between actual atom pos and that found from the use of F for all NN.
void TimeSeries::CalcLocalDefGradTensor(gsl_vector **NN_ti, gsl_vector **NN_tf, gsl_vector *w, int nNN, gsl_matrix *F)
{
	if (F==0){
		cout<<"Local Def Grad Tensor is not allocated."<<endl;
		exit(0);
	}
	gsl_matrix *D = gsl_matrix_calloc(3,3);
	gsl_matrix *A = gsl_matrix_calloc(3,3);
	//Calc D
	for (int i=0;i<nNN;i++){
		gsl_matrix *tempD=gsl_matrix_calloc(3,3);
		OuterProduct_gsl(NN_ti[i],NN_ti[i], tempD);
		double *w_i=gsl_vector_ptr(w,i);
		gsl_matrix_scale(tempD, *w_i);
		gsl_matrix_add(D,tempD);
		gsl_matrix_free(tempD);
	}
	//Calc A
	for (int i=0;i<nNN;i++){
		gsl_matrix *tempA=gsl_matrix_calloc(3,3);
		OuterProduct_gsl(NN_tf[i],NN_ti[i], tempA);
		double *w_i=gsl_vector_ptr(w,i);
		gsl_matrix_scale(tempA, *w_i);
		gsl_matrix_add(A,tempA);
		gsl_matrix_free(tempA);
	}
	//Invert D
	int s;
	gsl_permutation * p= gsl_permutation_alloc(3);
	gsl_linalg_LU_decomp(D,p,&s);
	gsl_matrix *DInv = gsl_matrix_calloc(3,3);
	gsl_linalg_LU_invert(D,p,DInv);
	// Calculate F  (AD^-1)
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A, DInv, 0.0, F);
	//Cleanup
	gsl_matrix_free(DInv);
	gsl_matrix_free(D);
	gsl_matrix_free(A);
	gsl_permutation_free(p);
}
void TimeSeries::CalcCGStrainTensor(gsl_matrix *F, gsl_matrix *A, string convention)
{
	if (F==0){
		cout<<"Error: Deformation Gradient Tensor is not allocated.\n";
		exit(0);
	}
	if (A==0){
		cout<<"Error: Left C-G Tensor is not allocated.\n";
		exit(0);
	}
	if (convention=="left"){
		gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, F, F, 0.0, A);
	}
	else if (convention=="right"){
		gsl_blas_dgemm(CblasTrans,CblasNoTrans, 1.0, F, F, 0.0, A);
	}
	else{
		cout<<"Convention: "<<convention<<" not supported in CalcCGStrainTensor."<<endl;
		exit(0);
	}
}

/* Decomposition of the Rotation matrix from the Deformation Gradent Tensor
 * 		Calculates Right and Left C-G Strain tensors
 * 		Determines their unit normal eigenvectors
 * 		Uses this to calculate the rotation matrix.
 * 		(see wikipedia: Finite Strain Theory for more info)
 */

void TimeSeries::DecompDefGradTensor(gsl_matrix *F, gsl_matrix *R, gsl_matrix *U, gsl_matrix *V)
{
	gsl_matrix *B=gsl_matrix_calloc(3,3), *C=gsl_matrix_calloc(3,3);
	CalcCGStrainTensor(F, C, "right");
	CalcCGStrainTensor(F,B, "left");
	//Alloc eigen vector computation
	/*
	gsl_vector_complex *lam_N=gsl_vector_complex_calloc(3);
	gsl_vector_complex *lam_n=gsl_vector_complex_calloc(3);
	gsl_matrix_complex *N=gsl_matrix_complex_calloc(3,3);
	gsl_matrix_complex *n=gsl_matrix_complex_calloc(3,3);
	gsl_eigen_nonsymmv_workspace *w1=gsl_eigen_nonsymmv_alloc(3);
	gsl_eigen_nonsymmv_workspace *w2=gsl_eigen_nonsymmv_alloc(3);
	*/
	gsl_vector *lam_N=gsl_vector_calloc(3);
	gsl_vector *lam_n=gsl_vector_calloc(3);
	gsl_matrix *N=gsl_matrix_calloc(3,3);
	gsl_matrix *n=gsl_matrix_calloc(3,3);
	gsl_eigen_symmv_workspace *w1=gsl_eigen_symmv_alloc(3);
	gsl_eigen_symmv_workspace *w2=gsl_eigen_symmv_alloc(3);
	//Calc eigenvectors of Right and Left C-G tensor
	//gsl_eigen_nonsymmv(C, lam_N, N, w1);
	//gsl_eigen_nonsymmv(B, lam_n, n, w2);
	gsl_eigen_symmv(C, lam_N, N, w1);
	gsl_eigen_symmv(B, lam_n, n, w2);
	gsl_eigen_symmv_sort(lam_N, N, GSL_EIGEN_SORT_VAL_DESC);
	gsl_eigen_symmv_sort(lam_n, n, GSL_EIGEN_SORT_VAL_DESC);
	/* Output the eigenvectors and eigenvalues
	OutputEigenVectors(outPath, name+"Eigen","EigenVals/Vectors of C", lam_N, N );
	OutputEigenVectors(outPath, name+"Eigen","EigenVals/Vectors of B", lam_n, n );
	OutputStrainTensor(&C, 1, outPath, name+"CGTens", "C Tensor");
	OutputStrainTensor(&B, 1, outPath, name+"CGTens", "B Tensor");
	*/
	//Cleanup
	//gsl_eigen_nonsymmv_free(w1);
	//gsl_eigen_nonsymmv_free(w2);
	gsl_eigen_symmv_free(w1);
	gsl_eigen_symmv_free(w2);
	//Alloc Temp matrices
	gsl_matrix *RTemp=gsl_matrix_calloc(3,3);
	gsl_matrix *UTemp=gsl_matrix_calloc(3,3);
	gsl_matrix *VTemp=gsl_matrix_calloc(3,3);
	gsl_matrix_set_zero(R);
	gsl_matrix_set_zero(U);
	gsl_matrix_set_zero(V);
	for (int i=0;i<3;i++){
		gsl_matrix_set_zero(RTemp);
		gsl_matrix_set_zero(UTemp);
		gsl_matrix_set_zero(VTemp);
		/*
		gsl_vector *Nreal=gsl_vector_calloc(3);
		gsl_vector *nreal=gsl_vector_calloc(3);
		for (int j=0;j<3;j++){
			//Check for imaginary part of determined eigenvectors
			gsl_complex *temp=gsl_matrix_complex_ptr(N, j, i);
			if(abs(GSL_IMAG(*temp))>PREC){
				cout<<"Nonreal eigenvectors found in DecompRotFromDefGradTensor.\n"<<endl;
				exit(0);
			}
			gsl_vector_set(Nreal, j, GSL_REAL(*temp));
			temp=gsl_matrix_complex_ptr(n, j, i);
			if(abs(GSL_IMAG(*temp))>PREC){
				cout<<"Nonreal eigenvectors found in DecompRotFromDefGradTensor.\n"<<endl;
				exit(0);
			}
			gsl_vector_set(nreal, j, GSL_REAL(*temp));
		}
		*/
		//Calc temporary outer product matrix
		gsl_vector_view n_col=gsl_matrix_column(n,i);
		gsl_vector_view N_col=gsl_matrix_column(N,i);
		OuterProduct_gsl(&n_col.vector, &N_col.vector, RTemp);
		OuterProduct_gsl(&n_col.vector, &n_col.vector, VTemp);
		OuterProduct_gsl(&N_col.vector, &N_col.vector, UTemp);
		//Calc Resulting Rotation matrix
		gsl_matrix_add(R, RTemp);
		//Calc Right stretch Tensor
		gsl_matrix_scale(UTemp, gsl_vector_get(lam_N, i));
		gsl_matrix_add(U, UTemp);
		//Calc Left Stretch Tensor
		gsl_matrix_scale(VTemp, gsl_vector_get(lam_n, i));
		gsl_matrix_add(V, VTemp);

		//Cleanup
		//gsl_vector_free(Nreal);
		//gsl_vector_free(nreal);
	}
	gsl_vector_free(lam_n);
	gsl_matrix_free(n);
	gsl_vector_free(lam_N);
	gsl_matrix_free(N);
	gsl_matrix_free(RTemp);
	gsl_matrix_free(UTemp);
	gsl_matrix_free(VTemp);
	gsl_matrix_free(B);
	gsl_matrix_free(C);
}
void TimeSeries::OutputEigenVectors(string path, string _name, string note, gsl_vector_complex *eVal, gsl_matrix_complex *eVect)
{
	string _file=path+"/"+_name+".evect";
	cout<< "\nOutputting Eigen Vectors to File: "<<_file<<endl;
	fstream file(_file.c_str(), fstream::out|fstream::app);

	file << note<<endl;
	file<< "EigenVal_1 EigenVect_1[1] EigenVect_1[2] EigenVect_1[3] ..."<<endl;;
	for (int i=0;i<3;i++){
		gsl_complex *com;
		com=gsl_vector_complex_ptr(eVal,i);
		file<<GSL_REAL(*com)<<" + i"<<GSL_IMAG(*com)<<"   ";
		for (int j=0;j<3;j++){
			com=gsl_matrix_complex_ptr(eVect,j,i);
			file<<GSL_REAL(*com)<<" + i"<<GSL_IMAG(*com)<<"   ";
		}
		file<<endl;
	}
	file.close();
}
void TimeSeries::OutputEigenVectors(string path, string _name, string note, gsl_vector *eVal, gsl_matrix *eVect)
{
	string _file=path+"/"+_name+".evect";
	cout<< "\nOutputting Eigen Vectors to File: "<<_file<<endl;
	fstream file(_file.c_str(), fstream::out|fstream::app);

	file << note<<endl;
	file<< "EigenVal_1 EigenVect_1[1] EigenVect_1[2] EigenVect_1[3] ..."<<endl;;
	for (int i=0;i<3;i++){
		file<<gsl_vector_get(eVal,i)<<" ";
		for (int j=0;j<3;j++){
			file<<gsl_matrix_get(eVect,i,j)<<" ";
		}
		file<<endl;
	}
	file.close();
}
void TimeSeries::ApplyHeavysideStepWeights(gsl_vector *w, gsl_vector **r, int n, double r_c)
{
	double D2=r_c*r_c;
	for (int i=0;i<n;i++){
		double d=0;
		gsl_blas_ddot(r[i], r[i], &d);
		double temp;
		if(d<D2){
			temp=1;
		}
		else if (d==D2){
			temp=0.5;
		}
		else{
			temp=0;
		}
		gsl_vector_set(w,i,temp);
	}
}

void TimeSeries::CalcLGStrainTensor(int s, AtomStrain *E)
{
	E->nA=F[s].nA;
	try{
		E->s=new gsl_matrix *[E->nA];
		E->index=new int[E->nA];
	}
	catch (exception &e){
	}
	for (int i=0;i<F[s].nA;i++){
		E->index[i]=F[s].index[i];
		//Calc E = 1/2(FT*F -I)
		E->s[i]=gsl_matrix_alloc(3,3);
		gsl_matrix_set_identity(E->s[i]);
		gsl_blas_dgemm(CblasTrans, CblasNoTrans, 0.5, F[s].s[i], F[s].s[i], -0.5, E->s[i]);
	}
}
void TimeSeries::CalcEAStrainTensor(int j, AtomStrain *e)
{
	e->nA=F[j].nA;
	try{
		e->s=new gsl_matrix *[e->nA];
		e->index=new int[e->nA];
	}
	catch (exception &f){
	}
	for (int i=0;i<F[j].nA;i++){
		e->index[i]=F[j].index[i];
		//Calc F*FT
		gsl_matrix *FFt=gsl_matrix_calloc(3,3);
		gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, F[j].s[i], F[j].s[i], 0, FFt);
		//Calc Inverse
		int s;
		gsl_permutation * p= gsl_permutation_alloc(3);
		gsl_linalg_LU_decomp(FFt,p,&s);
		gsl_matrix *FFtInv = gsl_matrix_calloc(3,3);
		gsl_linalg_LU_invert(FFt,p,FFtInv);
		//Calc e = 1/2(I- (F*FT)^-1)
		e->s[i]=gsl_matrix_alloc(3,3);
		gsl_matrix_set_identity(e->s[i]);
		gsl_matrix_sub(e->s[i],FFtInv);
		gsl_matrix_scale(e->s[i],0.5);
		//Cleanup
		gsl_matrix_free(FFt);
		gsl_matrix_free(FFtInv);
		gsl_permutation_free(p);
	}
}
void TimeSeries::CalcVelCorrFunc(string groupName, string tRName, double tScale)
{
	if (tScale<=0){
		cout<<endl<<"Error in CalcVelCorrFunc: Unphysical value of time scale given: "<<tScale<<endl;
		exit(0);
	}
	Range *tR=FindTimeRange(tRName);
	if (!tR){
		exit(0);
	}
	int del=tR->dx;
	if (del==0){
		del=GetSmallConstTimeStep(tR->x0, tR->xf, "vel");
	}
	if (del!=0){
		cout<<"\nCalculating Velocity Correlation Function ... ";
		int n=GetNumVelInRange(tR);
		AllocDistrib(&velCorrFn, "VelCorrFn", "Delta_t", "VelCorrelation", -(n-1)*del*tScale, (n-1)*del*tScale, del*tScale);
		int t0=tR->x0;
		int count=0;
		for (int i=0;i<n;i++){
			Timestep *time1=FindTimestep(t0);
			if (time1==0){
				cout<<"Defined timeRange must be exact, and timesteps must be regularly spaced.\n";
				exit(0);
			}
			time1->InputVel("all");
			velCorrFn->y[n-1]+=CalcAvgVelAutoCorr(groupName, time1, time1)/(double)n;
			for (int j=1;j<n-count;j++){
				int t2=t0+j*del;
				Timestep *time2=FindTimestep(t2);
				if (time2==0){
					cout<<"Defined timeRange must be exact, and timesteps must be regularly spaced.\n";
					exit(0);
				}
				time2->InputVel("all");
				velCorrFn->y[n-1+j]+=CalcAvgVelAutoCorr(groupName, time1, time2)/(double)(n-j);
				time2->Clean();
			}
			time1->Clean();
			t0+=del;
			count++;
		}
		//Make symmetrical distribution
		for(int i=0;i<n-1;i++){
			velCorrFn->y[i]=velCorrFn->y[2*(n-1)-i];
		}
		cout<<"done.\n";
	}
	else{
		cout<<"Error found velocity timestep of 0, make sure velocities are input. \n";
		exit(0);
	}
}
double TimeSeries::CalcAvgVelAutoCorr(string groupName, Timestep *v1, Timestep *v2)
{
	Group *g1=v1->FindGroup(groupName), *g2=v2->FindGroup(groupName);
	if (g1->nAtoms!=g2->nAtoms){
		cout<< "Cannot calculate Velocity Correlation because timeStep: "<<v1->step<<" number of atoms does not equal that of timeStep: " <<v2->step<<endl;
		exit(0);
	}
	//g1->CalcVelMags();
	//g2->CalcVelMags();
	int n=0;
	double sumVSq=0;
	for (int i=0;i<g1->nAtoms;i++){
		if (g1->atom[i]->index==g2->atom[i]->index){
			double vSq;
			gsl_blas_ddot(g1->atom[i]->vel,g2->atom[i]->vel, &vSq);
			sumVSq+=vSq;
			n++;
		}
		else{
			cout<<"Warning index of atoms not equal for CalcAvgVelAutoCorr: "<<i<<endl;
		}
	}
	sumVSq/=(double)n;
	return sumVSq;
}
void TimeSeries::AllocAvg(string groupName)
{
	try{
		if (avg){
			delete avg;
		}
		avg=new Group();
	}
	catch(exception &e){
		cout<<"Exception caught in TimeSeries::AllocAvg: "<<e.what()<<endl;
		exit(1);
	}
	avg->name=groupName+"Avg";
	avg->type="all";
}
//TODO Consider CalcAvgPos for group of atoms, not whole
void TimeSeries::CalcAvgPos(string groupName, int t0, int tf)
{
	cout<<"\nCalculating Average Position... ";
	AllocAvg(groupName);
	int n=0;
	for (int i=0;i<nTimeSteps;i++){
		if (step[i]->step>=t0&&step[i]->step<=tf&&step[i]->firstPos!=0){
			step[i]->InputPos("all");
			Group *g=step[i]->FindGroup(groupName);
			if (g->atom){
				if (!avg->atom){
					avg->CopyAtomInfo(g);
					avg->CopyAtomPos(g);
					try{
						if (avg->p){
							delete avg->p;
						}
						avg->p=new Position();
					}
					catch(exception &e){
					}
					if (g->p){
						for (int j=0;j<3;j++){
							for (int k=0;k<2;k++){
								avg->p->dim[j][k]=g->p->dim[j][k];
							}
							avg->p->sideL[j]=g->p->sideL[j];
						}
					}
				}
				else{
					if (avg->nAtoms!=g->nAtoms){
						cout<<"Cannot calculate average position, number of atoms are not equal,\n";
						exit(0);
					}
					//TODO Check atom pos for PBCs
					for (int i=0;i<g->nAtoms;i++){
						if (avg->atom[i]->index==g->atom[i]->index){
							gsl_vector_add(avg->atom[i]->pos, g->atom[i]->pos);
						}
						else{
							cout<<"Cannot calculate average position, atom indices are not equal.\n";
							exit(0);
						}
					}
					if (g->p){
						for (int j=0;j<3;j++){
							avg->p->dim[j][0]=getmin(avg->p->dim[j][0],g->p->dim[j][0]);
							avg->p->dim[j][1]=getmax(avg->p->dim[j][1],g->p->dim[j][1]);
							avg->p->sideL[j]=avg->p->dim[j][1]-avg->p->dim[j][0];
						}
					}
				}
				n++;
			}
			//clean up after calculation is finished
			step[i]->Clean();
		}
	}
	gsl_vector *norm=gsl_vector_alloc(3);
	gsl_vector_set_all(norm, 1.0/(double)n);
	avg->ScalePos(norm);
	gsl_vector_free(norm);
	cout<<"done.\n";
}
int TimeSeries::GetNumPosInRange(int t0, int tf)
{
	return GetNumPosInRange(t0,tf,0);
}
int TimeSeries::GetNumPosInRange(int t0, int tf, int dt)
{
	int count=0;
	for (int i=0;i<nTimeSteps;i++){
		if (step[i]->step>=t0&&step[i]->step<=tf){
			if(step[i]->firstPos!=0){
				count++;
				t0+=dt;
			}
		}
	}
	return count;
}
int TimeSeries::GetNumVelInRange(Range *tR)
{
	int count=0;
	int nextStep=tR->x0;
	for (int i=0;i<nTimeSteps;i++){
		if (step[i]->step>=nextStep && step[i]->step<=tR->xf){
			if (step[i]->firstVel!=0){
				count++;
				nextStep+=tR->dx;
			}
		}
	}
	return count;
}
int TimeSeries::GetMaxTimestep()
{
	int max=0;
	for (int i=0;i<nTimeSteps;i++){
		max=getmax(max,step[i]->step);
	}
	return max;
}
Timestep* TimeSeries::FindTimestep(int t)
{
	for (int i=0;i<nTimeSteps;i++){
		if (step[i]->step==t){
			return step[i];
		}
	}
	cout<<"Cannot Find Timestep "<<t<<" in sequence.\n";
	return 0;
}
void TimeSeries::CalcAvgAtomStrain(int j, gsl_matrix *A)
{
	//Average strain tensors of all atoms in F
	gsl_matrix_set_zero(A);
	if (F){
		for (int i=0;i<F[j].nA;i++){
			gsl_matrix_add(A,F[j].s[i]);
		}
		gsl_matrix_scale(A, 1/(double)F[j].nA);
	}
}
void TimeSeries::CleanAtomTensor(AtomStrain *A)
{
	try{
		if (A->s){
			for(int j=0;j<A->nA;j++){
				gsl_matrix_free(A->s[j]);
			}
			delete [] A->s;
		}
		if (A->index){
			delete [] A->index;
		}
	}
	catch(exception &e){

	}
}
void TimeSeries::OutputStrainTensor(AtomStrain *F, string path, string fileName, string comment)
{
	fileName=path+"/"+fileName;
	ofstream ofile;
//TODO check if strain tensor should be in scientific notation
	ofile.setf(ios::fixed, ios::floatfield);
	ofile.open(fileName.c_str(), ios_base::app);
	if (ofile.is_open()){
		ofile.precision(5);
		ofile<<comment<<endl;
		ofile<<"i F(0,0) F(0,1) F(0,2) F(1,0) F(1,1) F(1,2) F(2,0) F(2,1) F(2,2)\n";
		for (int i=0;i<F->nA;i++){
			ofile<<F->index[i]<<" ";
			for (int j=0;j<3;j++){
				for(int k=0;k<3;k++){
					ofile<<F->s[i]->data[j*F->s[i]->tda+k]<<" ";
				}
			}
			ofile <<endl;
		}
		ofile.close();
	}
	else{
		printf("\nCannot open file: %s\n",fileName.c_str());
		exit(0);
	}
	cout<<"Output tensor file: "<<fileName<<endl;
}
Range* TimeSeries::CreateNewTimeRange()
{
	Range **temp=0;
	if (tRange){
		if (nTimeRange>1){
			temp=new Range* [nTimeRange];
			for(int i=0;i<nTimeRange;i++){
				temp[i]=tRange[i];
			}
			delete [] tRange;
		}
		else{
			temp=new Range *;
			*temp=*tRange;
			delete tRange;
		}
	}
	nTimeRange++;
	if (nTimeRange>1){
		tRange=new Range* [nTimeRange];
		for (int i=0;i<nTimeRange-1;i++){
			tRange[i]=temp[i];
		}
		if (nTimeRange==2){
			delete temp;
		}
		else{
			delete [] temp;
		}
		tRange[nTimeRange-1]=new Range;
		return tRange[nTimeRange-1];
	}
	else{
		try{
			tRange = new Range *;
			*tRange=new Range;
		}
		catch (exception &e){
			cout<<e.what()<<endl;
		}
		return *tRange;
	}
}
Range* TimeSeries::FindTimeRange(string id)
{
	for (int i=0;i<nTimeRange;i++){
		//if (nTimeRange>1){
			if (tRange[i]->name==id){
				return tRange[i];
			}
		/*}
		else{
			if (*(tRange)->name==id){
				return *tRange;
			}
		}*/
	}
	cout<<"\nTimeRange '"<<id<<"' has not been defined.\n";
	return 0;
}

void TimeSeries::CalcNN_t(Distrib *dist, string tRName, string sRName, string g1Name, string g2Name)
{
	Range *tR=FindTimeRange(tRName);
	Range *sR=sRList.FindSpaceRange(sRName);
	int nTime=GetNumPosInRange(tR->x0, tR->xf, tR->dx);
	AllocDistrib(dist, "NeighDist(t)", "t", "d", 0, nTime, 1);
	int t0=tR->x0;
	int tf=tR->xf;
	int t=0;
	gsl_vector *d=gsl_vector_calloc(3);
	for (int i=0; i<nTimeSteps; i++){
		if (step[i]->step>=t0&& step[i]->step<=tf&&step[i]->firstPos!=0){
			double dTot=0, count=0;
			step[i]->InputPos("all");
			Group *g1=step[i]->FindGroup(g1Name);
			Group *g2=step[i]->FindGroup(g2Name);
			//g1->BuildNeighsInRange(g2, sR->x0, sR->xf);
			g1->BuildNeighsInRange_fast(g2, sR->x0, sR->xf);
			//double lowBound=sR->x0*sR->x0, highBound=sR->xf*sR->xf;
			double L_2[3];
			for(int j=0;j<3;j++){
				L_2[j]=g1->p->sideL[j]/2.0;
			}
			for (int j=0;j<g1->nAtoms;j++){
				for (int k=0; k<g1->atom[j]->n_NN;k++){
					CalcDistPBC_gsl(g1->atom[j]->NN[k]->pos, g1->atom[j]->pos, d, &g1->p->boundary[0], &L_2[0], &g1->p->sideL[0]);
					double temp;
					gsl_blas_ddot(d, d, &temp);
					dTot+=sqrt(temp);
					count++;
					/*
					double temp=0;
					gsl_vector_memcpy(d, g1->atom[j]->NN[k]->pos);
					gsl_vector_sub(d, g1->atom[j]->pos);
					gsl_blas_ddot(d, d, &temp);
					if (temp>=lowBound&&temp<=highBound){
						dTot+=sqrt(temp);
						count++;
					}
					*/
				}
			}
			dist->x[t]=step[i]->step;
			dist->y[t]=dTot/count;
			t++;
			step[i]->Clean();
			t0=step[i]->step+tR->dx;
		}
	}
	gsl_vector_free(d);
}
void TimeSeries::MovePosToNewTime(int tOld, string oldName, int tNew, string newName)
{
	Timestep *oldTime=FindTimestep(tOld);
	Timestep *newTime=FindTimestep(tNew);
	Group *oldPos=0;
	if (!oldTime){
		exit(0);
	}
	else{
		oldPos=oldTime->FindPos(oldName);
		if (!oldPos){
			cout<<"Cannot find position with name: "<<oldName<<endl;
		}
	}
	if (!newTime){
		cout<<"Creating new Timestep: "<<tNew<<endl;
		Timestep *zero=0;
		newTime=firstStep;
		InsertStepInSeq(&zero,&newTime,tNew);
	}
	else{
		Group *temp=newTime->FindPos(newName);
		if(temp!=0){
			cout<<"Position name: "<<newName<<" already exists... exiting.\n";
			exit(0);
		}
		else{
			cout<<"Moving position to new timestep: "<<tNew<<endl;
		}
	}
	newTime->step=tNew;
	//add position to list of datafiles for timestep
	Group *newPos=new Group;
	if (newTime->curPos){
		newTime->curPos->next=newPos;
	}
	newTime->curPos=newPos;
	newTime->nPosData++;
	if (!newTime->firstPos){
		newTime->firstPos=newTime->curPos;
	}
	//
	newTime->curPos->AllocPos();
	newTime->curPos->name=newName;
	newTime->curPos->type="pos";
	newTime->curPos->p->file=oldPos->p->file;
	newTime->curPos->p->dataPos=oldPos->p->dataPos;
	newTime->curPos->p->fileType=oldPos->p->fileType;
	for (int i=0;i<3;i++){
		newTime->curPos->p->boundary[i]=oldPos->p->boundary[i];
	}
	newTime->CreatePosArray();
	oldTime->DeletePosData(oldName);
	//TODO write a check which then looks at timesteps and determines if some should be deleted.
	// They would only be kept if other information still exists (i.e. pos or vel)
	//Create indexed array from new linked list
	CreateSeriesArray();
}
void TimeSeries::DefineUserGroup(string tRName, string groupName, string groupComm)
{
	Range *tR=FindTimeRange(tRName);
	int t0=tR->x0;
	bool foundTimestep=false;
	for (int i=0;i<nTimeSteps;i++){
		if (step[i]->step >= t0 && step[i]->step<=tR->xf){
			foundTimestep=true;
			step[i]->AddUserGroup(groupName, groupComm);
			t0+=tR->dx;
		}
	}
	if (!foundTimestep){
		cout<<"Warning: Could not define group: "<<groupName<<", no timesteps found in range.\n";
	}
}
void TimeSeries::CalcRadius(string tRName, string groupName, string path, string fileName)
{
	Range *tR=FindTimeRange(tRName);
	int t0=tR->x0;
	//cout<<"Radius Info of "<<groupName<<" atoms:\n"<<"time min max avg var Rg\n";

	Distrib minR, maxR, avgR, varR, Rg;
	AllocDistrib(&minR, "minR", "time", "minR", tR->x0,tR->xf, tR->dx);
	AllocDistrib(&maxR, "maxR", "time", "maxR", tR->x0,tR->xf, tR->dx);
	AllocDistrib(&avgR, "avgR", "time", "avgR", tR->x0,tR->xf, tR->dx);
	AllocDistrib(&varR, "varR", "time", "varR", tR->x0,tR->xf, tR->dx);
	AllocDistrib(&Rg, "Rg", "time", "Rg", tR->x0,tR->xf, tR->dx);
	for (int i=0;i<nTimeSteps;i++){
		if (step[i]->step >= t0 && step[i]->step<=tR->xf){
			cout<<"Calculating radius for step: "<<step[i]->step<<endl;
			step[i]->InputPos("all");
			Group *g1=step[i]->FindGroup(groupName);
			step[i]->all.CalcCenter();
			g1->CalcPosMags(step[i]->all.p->center);
			double avg=g1->GetAvgR();
			minR.x[i]=step[i]->step;

			minR.y[i]=g1->GetMinR();
			maxR.y[i]=g1->GetMaxR();
			avgR.y[i]=avg;
			varR.y[i]=g1->GetVarR(avg);
			Rg.y[i]=g1->CalcRadiusOfGyration(step[i]->all.p->center);
			//cout<<step[i]->step<<" "<<g1->GetMinR()<<" "<<g1->GetMaxR()<<" "<<avg<<" "<<g1->GetVarR(avg)<<" "<<g1->CalcRadiusOfGyration(step[i]->all.p->center)<<endl;
			step[i]->Clean();
			t0=step[i]->step+tR->dx;
		}
	}

	Distrib **rDs=new Distrib *[5];
	rDs[0]=&minR; rDs[1]=&maxR; rDs[2]=&avgR; rDs[3]=&varR;; rDs[4]=&Rg;
	OutputMultDistribs(rDs, 5, path, fileName, false, false, true, "Radius Info of "+groupName+" atoms");
	delete [] rDs;

	CleanDistrib(&minR);
	CleanDistrib(&maxR);
	CleanDistrib(&avgR);
	CleanDistrib(&varR);
	CleanDistrib(&Rg);

}
