/*
 * Group.cpp
 *
 *  Created on: Nov 18, 2010
 *      Author: ken
 */

#include "Group.h"

Group::Group() {
	name="";
	type="";
	nAtoms=0;
	atom=0;
	p=0;
	v=0;
	next=0;
	neighR0=0;
	neighRc=0;
	neighG=0;
	nSubDomX=0;
	nSubDomY=0;
	nSubDomZ=0;
	subDomain=0;
	genCommand="";
}

Group::~Group() {
	if (type!="user"){
		delete p;
		delete v;
	}
	if (subDomain){
		for(int i=0;i<nSubDomX;i++){
			for (int j=0; j<nSubDomY; j++){
				for (int k=0;k<nSubDomZ;k++){
					delete [] subDomain[i][j][k].atom;
					subDomain[i][j][k].atom=0;
					subDomain[i][j][k].nAtoms=0;
				}
				delete [] subDomain[i][j];
			}
			delete [] subDomain[i];
		}
		delete [] subDomain;
		subDomain=0;
		nSubDomX=0;
		nSubDomY=0;
		nSubDomZ=0;
	}
	if (atom){
		for (int i=0;i<nAtoms;i++){
			if (type=="all"){
				if (atom[i]->pos){
					gsl_vector_free(atom[i]->pos);
				}
				if (atom[i]->vel){
					gsl_vector_free(atom[i]->vel);
				}
				if (atom[i]->NN){
					delete [] atom[i]->NN;
				}
				delete atom[i];
			}
			if (type=="pos"){
				if (atom[i]->pos){
					gsl_vector_free(atom[i]->pos);
				}
				delete atom[i];
			}
			if (type=="vel"){
				if (atom[i]->vel){
					gsl_vector_free(atom[i]->vel);
				}
				delete atom[i];
			}
		}
		delete [] atom;
		atom=0;
	}

}
void Group::AllocAtoms(){
	try{
		if (atom){
			for (int i=0;i<nAtoms;i++){
				if (atom[i]){
					if (type=="all"||type=="pos"){
						if (atom[i]->pos){
							gsl_vector_free(atom[i]->pos);
						}
					}
					if (type=="all"||type=="vel"){
						if (atom[i]->vel){
							gsl_vector_free(atom[i]->vel);
						}
					}
					delete atom[i];
				}
			}
			delete [] atom;
		}
		//only "pos", "vel" and "all" type groups carry the information, other groups are just pointers to the all group
		atom=new Atom *[nAtoms];
		for (int i=0;i<nAtoms;i++){
			atom[i]=0;
			if (type=="pos"||type=="vel"||type=="all"){
				atom[i]=new Atom;
				ZeroAtom(atom[i]);
			}
		}
	}
	catch(exception &e){
		cout<<"Exception caught in Group::AllocAtoms: "<<e.what()<<endl;
		exit(1);
	}
}
void Group::AllocPos(){
	try{
		if (p){
			delete p;
		}
		p=new Position;
	}
	catch(exception &e){
		cout<<"Exception caught in Group::AllocPos: "<<e.what()<<endl;
		exit(1);
	}
}
void Group::AllocVel(){
	try{
		if (v){
			delete v;
		}
		v=new Velocity;
	}
	catch(exception &e){
		cout<<"Exception caught in Group::AllocVel: "<<e.what()<<endl;
		exit(1);
	}
}
void Group::AllocAtomPos(){
	try{
		if (!atom){
			cout<<"Cannot Allocate Atom Positions, Atoms are not allocated.\n";
			exit(0);
		}
		for (int i=0;i<nAtoms;i++){
			if (atom[i]){
				if (atom[i]->pos){
					gsl_vector_free(atom[i]->pos);
				}
				//TODO Read up on and implement error handling with gsl handler
				atom[i]->pos=gsl_vector_calloc(3);
			}
			else{
				cout<<"Cannot Allocate Atom Positions, Atom "<<i<<" is not allocated.\n";
				exit(0);
			}
		}
	}
	catch(exception &e){
		cout<<"Exception caught in Group::AllocAtomPos: "<<e.what()<<endl;
		exit(1);
	}
}
void Group::AllocAtomVel(){
	try{
		if (!atom){
			cout<<"Cannot Allocate Atom Velocities, Atoms are not allocated.\n";
			exit(0);
		}
		for (int i=0;i<nAtoms;i++){
			if (atom[i]){
				if (atom[i]->vel){
					gsl_vector_free(atom[i]->vel);
				}
				//TODO Read up on and implement error handling with gsl handler
				atom[i]->vel=gsl_vector_calloc(3);
			}
			else{
				cout<<"Cannot Allocate Atom Positions, Atom "<<i<<" is not allocated.\n";
				exit(0);
			}
		}
	}
	catch(exception &e){
		cout<<"Exception caught in Group::AllocAtomVel: "<<e.what()<<endl;
		exit(1);
	}
}
void Group::SortByIndex()
{
	//Copy important info to arrays and sort arrays
	double *ind, *place;
	try{
		ind=new double [nAtoms];
		place=new double [nAtoms];
	}
	catch(exception &e){
		cout<<"Exception in SortByIndex: "<<e.what()<<endl;
		exit(1);
	}
	for (int i=0;i<nAtoms;i++){
		ind[i]=(double)atom[i]->index;
		place[i]=(double)i;
	}
	QuickSortImproved(ind, place, 0, nAtoms-1);

	//Make sorted atom array
	Atom **tempAtom;
	try{
		tempAtom=new Atom* [nAtoms];
	}
	catch(exception &e){
		cout<<"Exception in SortByIndex: "<<e.what();
		exit (1);
	}
	for (int i=0;i<nAtoms;i++){
		tempAtom[i]=atom[(int)place[i]];
	}
	//Set atom pointer to sorted list of atoms
	delete [] atom;
	atom=tempAtom;

	delete [] ind;
	delete [] place;
}
void Group::CopyAtomInfo(Group *a)
{
	if (nAtoms==0){
		nAtoms=a->nAtoms;
		AllocAtoms();
	}
	for (int i=0;i<nAtoms;i++){
		atom[i]->index=a->atom[i]->index;
		atom[i]->type=a->atom[i]->type;
	}
}
void Group::CopyAtomPos(Group *a)
{
	if (nAtoms==0){
		nAtoms=a->nAtoms;
		AllocAtoms();
	}
	AllocAtomPos();
	for (int i=0;i<nAtoms;i++){
		gsl_vector_memcpy(atom[i]->pos, a->atom[i]->pos);
		atom[i]->r=a->atom[i]->r;
	}
}
void Group::CopyAtomVel(Group *a)
{
	if (nAtoms==0){
		nAtoms=a->nAtoms;
		AllocAtoms();
	}
	AllocAtomVel();
	for (int i=0;i<nAtoms;i++){
		gsl_vector_memcpy(atom[i]->vel, a->atom[i]->vel);
		atom[i]->v=a->atom[i]->v;
	}
}
void Group::CalcVelDistrib()
{
	double con=1.0/(double)nAtoms;
	double avgVel=0;
	CalcVelMags();
	if (!v->pVel){
		AllocDistrib(&v->pVel, "pVel", "Velocity", "pVel", 0, GetMaxVel(), VPREC);
	}
	for (int i=0;i<nAtoms;i++){
		avgVel+=atom[i]->v;
		v->pVel->y[(int)(atom[i]->v/v->pVel->step)]+=con;
	}
	avgVel/=(double)nAtoms;
	cout<<"Avg. atomic velocity of group "<<name<< ": "<<avgVel<<endl;
}
void Group::CalcVelMags()
{
	for (int i=0; i<nAtoms; i++){
		if (atom[i]->v==0){
			gsl_blas_ddot(atom[i]->vel, atom[i]->vel, &atom[i]->v);
			atom[i]->v=sqrt(atom[i]->v);
		}
	}
}
void Group::CalcPosMags(gsl_vector *zero)
{
	for (int i=0; i<nAtoms; i++){
		gsl_vector *temp=gsl_vector_alloc(3);
		gsl_vector_memcpy(temp,atom[i]->pos);
		gsl_vector_sub(temp, zero);
		gsl_blas_ddot(temp, temp, &atom[i]->r);
		atom[i]->r=sqrt(atom[i]->r);
		gsl_vector_free(temp);
	}
}
double Group::GetMaxVel()
{
	double maxV=0;
	CalcVelMags();
	for (int i=0;i<nAtoms;i++){
		maxV=getmax(atom[i]->v,maxV);
	}
	return maxV;
}
double Group::GetMaxR()
{
	double maxR=0;
	for (int i=0;i<nAtoms;i++){
		if (atom[i]->r!=-1){
			maxR=getmax(atom[i]->r,maxR);
		}
		else{
			cout<<"Error: cannot get max of radial pos in group: "<<name<<", radial pos is not yet defined.\n";
			exit(0);
		}

	}
	return maxR;
}
double Group::GetMinR()
{
	double minR;
	for (int i=0;i<nAtoms;i++){
		if (atom[i]->r!=-1){
			if (i==0){
				minR=atom[i]->r;
			}
			else{
				minR=getmin(atom[i]->r,minR);
			}
		}
		else{
			cout<<"Error: cannot get max of radial pos in group: "<<name<<", radial pos is not yet defined.\n";
			exit(0);
		}
	}
	return minR;
}
double Group::GetAvgR()
{
	double avgR=0;
	for (int i=0;i<nAtoms;i++){
		if (atom[i]->r!=-1){
			avgR+=atom[i]->r;
		}
		else{
			cout<<"Error: cannot calc avg of radial pos in group: "<<name<<", radial pos is not yet defined.\n";
			exit(0);
		}
	}
	return avgR/(double)nAtoms;
}
double Group::GetVarR(double avgR)
{
	double var=0;
	for (int i=0;i<nAtoms;i++){
		if (atom[i]->r!=-1){
			double temp=(atom[i]->r-avgR);
			var+=temp*temp;
		}
		else{
			cout<<"Error: cannot calculate variance of radial pos in group: "<<name<<", radial pos is not yet defined.\n";
			exit(0);
		}
	}
	var/=(double)nAtoms;
	return var;
}
void Group::ScalePos(gsl_vector *scale)
{
	for (int i=0;i<nAtoms;i++){
		gsl_vector_mul(atom[i]->pos, scale);
	}
	if (p){
		if (p->center){
			gsl_vector_mul(p->center, scale);
		}
	}
}
void Group::TransPos(gsl_vector *dr)
{
	for (int i=0;i<nAtoms;i++){
		gsl_vector_add(atom[i]->pos,dr);
	}
	if (p->center){
		gsl_vector_add(p->center,dr);
	}
}
void Group::RotatePos(gsl_matrix *R, gsl_vector *pivot)
{
	for (int i=0;i<nAtoms;i++){
		gsl_vector *temp=gsl_vector_alloc(3);
		gsl_vector_memcpy(temp, atom[i]->pos);
		gsl_vector_sub(temp, pivot);
		//TODO Debug here
		gsl_blas_dgemv (CblasNoTrans, 1.0, R, temp, 0, temp);
		gsl_vector_add(temp, pivot);
		gsl_vector_memcpy(atom[i]->pos, temp);
	}
}
void Group::RemapPos(int dim, double min, double max)
{
	double L=max-min;
	for (int i=0;i<nAtoms;i++){
		if (atom[i]->pos->data[dim]<min){
			atom[i]->pos->data[dim]+=L;
		}
		if (atom[i]->pos->data[dim]>=max){
			atom[i]->pos->data[dim]-=L;
		}
	}
}
void Group::CalcCenter()
{
	if (!p->center){
		p->center=gsl_vector_calloc(3);
	}
	else {
		gsl_vector_set_zero(p->center);
	}
	for (int i=0;i<nAtoms;i++){
		gsl_vector_add(p->center, atom[i]->pos);
	}
	gsl_vector_scale(p->center,1/(double) nAtoms);
}
int Group::GetTotNumNN()
{
	int tot_nNN=0;
	for (int i=0;i<nAtoms;i++){
		tot_nNN+=atom[i]->n_NN;
	}
	return tot_nNN;
}
void Group::AllocNN()
{
	for (int i=0;i<nAtoms;i++){
		AllocAtomNN(atom[i]);
	}
}
void Group::CreateRadialDistrib(gsl_vector *center, double prec)
{
	if (!p){
		p=new Position;
	}
	//TODO Consider how to make work without recalculating everything
	AllocDistrib(&p->rad, "RadialDensity", "RadialPosition","N", 0,p->GetMaxD(), prec);
	CalcPosMags(center);
	for (int i=0;i<nAtoms;i++){
		p->rad->y[(int)(atom[i]->r/p->rad->step)]++;
	}
}
void Group::OutputNumNNLammpstrjFile(string path, string simName, int t)
{
	string _file=path+"/"+simName+name+"NN.lammpstrj";
	cout<< "\nOutputting NN Particle Position to File: "<<_file<<endl;
	fstream file(_file.c_str(), fstream::out|fstream::app);

	file << "ITEM: TIMESTEP\n"<<t<<"\n" << "ITEM: NUMBER OF ATOMS\n" << nAtoms << "\n";
	file << "ITEM: BOX BOUNDS\n" ;
	for (int i=0;i<3;i++){
		for(int j=0;j<2;j++){
			file<< p->dim[i][j] << " ";
		}
		p->sideL[i]=p->dim[i][1]-p->dim[i][0];
		file<<endl;
	}
	file << "ITEM: ATOMS\n";
	for (int i=0;i<nAtoms;i++){
		file.precision(0);
		file << atom[i]->index << " " << atom[i]->n_NN << " ";
		file.precision(5);
		for (int k=0;k<3;k++){
			file<<(gsl_vector_get(atom[i]->pos,k)-p->dim[k][0])/p->sideL[k]<< " " ;
		}
		file<<endl;
	}
	file.close();
}
void Group::OutputNNVectorLammpstrjFile(string path, string simName, int t)
{
	string _file=path+"/"+simName+name+"NNVect.lammpstrj";
	cout<< "\nOutputting NN Vectors to File: "<<_file<<endl;
	fstream file(_file.c_str(), fstream::out|fstream::app);

	file << "ITEM: TIMESTEP\n"<<t<<"\n" << "ITEM: NUMBER OF ATOMS\n" << GetTotNumNN() << "\n";
	file << "ITEM: BOX BOUNDS\n0 1\n0 1\n0 1\n" ;
	file << "ITEM: ATOMS\n";
	int j=1;
	for (int i=0;i<nAtoms;i++){
		for (int k=0;k<atom[i]->n_NN;k++){
			file.precision(0);
			file << j << " " << atom[i]->n_NN<< " ";
			file.precision(5);
			gsl_vector *d=gsl_vector_alloc(3);
			gsl_vector_memcpy(d,atom[i]->NN[k]->pos);
			gsl_vector_sub(d,atom[i]->pos);
			for (int dim=0;dim<3;dim++){
				file<< gsl_vector_get(d,dim)<< " " ;
			}
			file<<endl;
			j++;

		}
	}
	file.close();
}
void Group::CalcRDF(double prec)
{
	AllocDistrib(&p->rdf, "RDF", "d", "N", 0, p->GetMaxD(), prec);
	gsl_vector *temp=gsl_vector_alloc(3);
	double minSq=p->rdf->x[0]*p->rdf->x[0];
	double maxSq=p->rdf->x[p->rdf->n-1]*p->rdf->x[p->rdf->n-1];
	for (int i=0;i<nAtoms;i++){
		for (int j=i+1;j<nAtoms;j++){
			double d=0;
			gsl_vector_memcpy(temp, atom[i]->pos);
			gsl_vector_sub(temp,atom[j]->pos);
			gsl_blas_ddot(temp, temp, &d);
			if (d>minSq&&d<maxSq){
				p->rdf->y[(int)(sqrt(d)/p->rdf->step)]++;
			}
		}
	}
}
void Group::CalcRDF(Distrib *rdf)
{
	//gsl_vector *temp=gsl_vector_alloc(3);
	double minSq=rdf->x[0]*rdf->x[0];
	double maxSq=rdf->x[rdf->n-1];//+rdf->step/2.0;
	maxSq*=maxSq;
	for (int i=0;i<nAtoms;i++){
		for (int j=0;j<atom[i]->n_NN;j++){
			double d=0;
			gsl_vector *dist=gsl_vector_calloc(3);
			gsl_vector_memcpy(dist, atom[i]->NN[j]->pos);
			gsl_vector_sub(dist, atom[i]->pos);
			gsl_blas_ddot(dist, dist, &d);
			if (d>minSq&&d<maxSq){
				rdf->y[(int)((sqrt(d)-rdf->minX)/rdf->step)]++;
			}
			gsl_vector_free(dist);
		}
	}
	//gsl_vector_free(temp);
}
void Group::OutputVolChangeMap(double *dV, string path, string _name, int timeStep, int index)
{
	//TODO make output in PDB format (formatted so have to use setw() etc...
	double min, max, step;
	min=dV[0];
	max=dV[0];
	double sum=0;
	double sumSq=0;
	for (int i=0;i<nAtoms;i++ ){
		sum+=dV[i];
		sumSq+=dV[i]*dV[i];
		min=getmin(min, dV[i]);
		max=getmax(max, dV[i]);
	}
	double avg=sum/nAtoms;
	double avgSq=sumSq/nAtoms;
	double var= avgSq-avg*avg;
	double stdDev=sqrt(var);
	step=stdDev;
	string _file=path+"/"+_name;
	cout<< "\nOutputting Vol Change Map of timestep "<<timeStep<<" to File: "<<_file<<endl;
	cout<<"Volume Change Range: "<<min <<" "<<max<<" Map Step: "<<step<<endl;
	cout<<"--- Avg: "<< avg<<"StdDev: "<<stdDev<<endl;
	fstream file(_file.c_str(), fstream::out|fstream::app);

	//Align output to right of field and use fixed decimal precision
	file << right<<fixed;
	file <<"MODEL "<<setw(8)<<index%10000<<endl;
	for (int i=0;i<nAtoms;i++){
		file<<"ATOM  "<<setw(5)<<atom[i]->index%100000<<setw(4)<<'X'<<" "<<"VAL"<<" A"<<setw(4)<<"1"<<" "<<"   ";
		file.precision(3);
		for (int k=0;k<3;k++){
			file<<setw(8)<<gsl_vector_get(atom[i]->pos,k);
		}
		file << setw(6)<<setprecision(2)<< (dV[i]-avg)/stdDev ;
		file<<setw(6)<< "0.0"<<setw(6)<<" "<<setw(4)<<left<<"A1"<<right<<setw(2)<<"X"<<setw(2)<<"0";
		file<<endl;
	}
	file<<"ENDMDL"<<endl;
	file.close();
}
Atom* Group::FindAtom(int index)
{
	//First check members of the group
	for (int i=0;i<nAtoms;i++){
		if (atom[i]->index==index){
			return atom[i];
		}
	}
	//Then check members of Nearest Neighbors which might be outside group
	for (int i=0;i<nAtoms;i++){
		if (atom[i]->NN!=0){
			for (int j=0;j<atom[i]->n_NN;j++){
				if (atom[i]->NN[j]->index==index){
					return atom[i]->NN[j];
				}
			}
		}
	}
	cout<<"Cannot find atom: "<<index<<endl;
	exit(0);
}
void Group::OutputStrainTensor(gsl_matrix **F, int n, string path, string fileName, string comment)
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
		for (int i=0;i<n;i++){
			ofile<<atom[i]->index<<" ";
			for (int j=0;j<3;j++){
				for(int k=0;k<3;k++){
					ofile<<F[i]->data[j*F[i]->tda+k]<<" ";
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
void Group::BuildNeighsInRange(Group *g2, double R0, double Rc)
{
	//Check if calculation is between atoms in the same group
	bool sameGroup=(name==g2->name);
	if (!sameGroup){
		//If they are not the same group, double check that the sim box is the same.
		if (!p->SameBoxInfo(g2->p)){
			cout<<"Cannot Calc NN between group "<<this->name<<" and group "<<g2->name<<endl;
			cout<<"Not the same simulation box.\n";
			exit(0);
		}
	}
	double RcSq=Rc*Rc, R0Sq=R0*R0;
	double innerBound[3][2], outerBound[3][2];
	gsl_vector *temp=gsl_vector_alloc(3);
	int tempn_NN=0;
	Atom **tempNN=new Atom* [g2->nAtoms];
	//TODO consider removing zeroing and just copying nPointers
	//Create temporary NN list which is then copied over to each atom
	for (int i=0;i<g2->nAtoms;i++){
		tempNN[i]=0;
	}
	for (int i=0;i<nAtoms;i++){
		for (int j=0; j<tempn_NN;j++){
			tempNN[j]=0;
		}
		tempn_NN=0;
		// outerBound is cube inscribing a sphere of radius Rc
		// innerBound is cube inscribed in sphere of radius R0
		for (int dim=0;dim<3;dim++){
			double x=gsl_vector_get(atom[i]->pos, dim);
			outerBound[dim][0]=x-Rc;
			outerBound[dim][1]=x+Rc;
			innerBound[dim][0]=x-(R0/sqrt(3.0));
			innerBound[dim][1]=x+(R0/sqrt(3.0));
		}
		int start=0;

		if (sameGroup){
			//Check NN lists of already calculated atoms
			for (int j=0;j<i;j++){
				for (int k=0;k<atom[j]->n_NN;k++){
					if (atom[j]->NN[k]==atom[i]){
						tempNN[tempn_NN]=atom[j];
						tempn_NN++;
					}
				}
			}
			start=i;
		}

		//Find new NN
		for (int j=start; j<g2->nAtoms;j++){
			bool inBounds=true;
			// To diminish unnecessary distance calcs, create box around each atom and use inequalities.
			int nInner=0;
			double p2[3];
			for (int dim=0;dim<3&&inBounds;dim++){
				p2[dim]=gsl_vector_get(g2->atom[j]->pos,dim);

				if (p2[dim]<outerBound[dim][0]){
					if (p->boundary[dim]=='p'){
						if (outerBound[dim][1]>p->dim[dim][1]){
							inBounds=((p2[dim]+=p->sideL[dim])<outerBound[dim][1]);
						}
					}
					else{
						inBounds=false;
					}
				}
				else if (p2[dim]>outerBound[dim][1]){
					if (p->boundary[dim]=='p'){
						if (outerBound[dim][0]<p->dim[dim][0]){
							inBounds=((p2[dim]-=p->sideL[dim])>outerBound[dim][0]);
						}
					}
					else{
						inBounds=false;
					}
				}
				else{
					inBounds=true;
				}

				if (inBounds){
					if (p2[dim]>innerBound[dim][0]&&p2[dim]<innerBound[dim][1]){
						nInner++;
					}
				}
			}
			if (inBounds &&nInner==3){
				inBounds=false;
			}

			if (inBounds){
				// Positions already corrected for PBCs
				double dSq=0;

				gsl_vector_view p_2=gsl_vector_view_array(&p2[0],3);
				gsl_vector_sub(&p_2.vector, atom[i]->pos);
				gsl_blas_ddot(&p_2.vector, &p_2.vector, &dSq);

				if (dSq<=RcSq&&dSq>=R0Sq){
					tempNN[tempn_NN]=g2->atom[j];
					tempn_NN++;
				}
			}
		}
		//Copy NN list to atoms
		atom[i]->n_NN=tempn_NN;
		AllocAtomNN(atom[i]);
		for(int j=0;j<tempn_NN;j++){
			atom[i]->NN[j]=tempNN[j];
		}
	}
	gsl_vector_free(temp);
	delete [] tempNN;
}

void Group::BuildNeighsInRange_fast(Group *g2, double R0, double Rc)
{
	//Setup so SubDomain size is fraction of boxSize
	int n_L[3];
	for (int i=0;i<3;i++){
		n_L[i]=(int) (p->sideL[i]/Rc );
	}

	//Check if calculation is between atoms in the same group
	bool sameGroup=(name==g2->name);
	if (!sameGroup){
		//If they are not the same group, double check that the sim box is the same.
		if (!p->SameBoxInfo(g2->p)){
			cout<<"Cannot Calc NN between group "<<this->name<<" and group "<<g2->name<<endl;
			cout<<"Not the same simulation box.\n";
			exit(0);
		}
		g2->CreateSubDomains(&n_L[0]);
	}
	CreateSubDomains(&n_L[0]);

	double RcSq=Rc*Rc, R0Sq=R0*R0;
	double innerBound[3][2], outerBound[3][2];
	double innerL_2=R0/sqrt(3.0);
	int tempn_NN=0;
	Atom **tempNN=new Atom* [g2->nAtoms];
	//TODO consider removing zeroing and just copying nPointers
	//Create temporary NN list which is then copied over to each atom
	for (int i=0;i<g2->nAtoms;i++){
		tempNN[i]=0;
	}
	gsl_vector *d=gsl_vector_alloc(3);
	for (int i=0;i<nSubDomX;i++){
		for (int j=0;j<nSubDomY;j++){
			for (int k=0;k<nSubDomZ;k++){
				for (int a=0;a<subDomain[i][j][k].nAtoms;a++){
					//Zero tempNN List
					for (int x=0; x<tempn_NN;x++){
						tempNN[x]=0;
					}
					tempn_NN=0;
					// Define two limiting spheres for boundaries of what is considered a neighbor
					for (int dim=0;dim<3;dim++){
						double x=subDomain[i][j][k].atom[a]->pos->data[dim];
						outerBound[dim][0]=x-Rc;
						outerBound[dim][1]=x+Rc;
						innerBound[dim][0]=x-innerL_2;
						innerBound[dim][1]=x+innerL_2;
					}
					//Check all neighboring subDomains
					for (int l=i-1;l<=i+1;l++){
						for(int m=j-1;m<=j+1;m++){
							for (int n=k-1;n<=k+1;n++){
								bool outOfBox=false;
								// TODO cleanup the following lines of code since they are repetative
								int l2=l,m2=m,n2=n;
								if (p->boundary[0]=='p'){
									if (l==-1){
										l2+=nSubDomX;
									}
									if (l==nSubDomX){
										l2-=nSubDomX;
									}
								}
								else if (p->boundary[0]=='m'){
									if (l==-1||l==nSubDomX){
										outOfBox=true;
									}
								}
								else{
									cout<<"Error; boundary conditions not supported.\n";
									exit(0);
								}
								if (p->boundary[1]=='p'){
									if (m==-1){
										m2+=nSubDomY;
									}
									if (m==nSubDomY){
										m2-=nSubDomY;
									}
								}
								else if (p->boundary[1]=='m'){
									if (m==-1||m==nSubDomX){
										outOfBox=true;
									}
								}
								else{
									cout<<"Error; boundary conditions not supported.\n";
									exit(0);
								}
								if (p->boundary[2]=='p'){
									if (n==-1){
										n2+=nSubDomZ;
									}
									if (n==nSubDomZ){
										n2-=nSubDomZ;
									}
								}
								else if (p->boundary[2]=='m'){
									if (n==-1||n==nSubDomX){
										outOfBox=true;
									}
								}
								else{
									cout<<"Error; boundary conditions not supported.\n";
									exit(0);
								}
								if (!outOfBox){
									for (int b=0;b<g2->subDomain[l2][m2][n2].nAtoms;b++){
										bool inBounds=true;
										// To diminish unnecessary distance calcs, create box around each atom and use inequalities.
										int nInner=0;
										double p2[3];
										for (int dim=0;dim<3&&inBounds;dim++){
											p2[dim]=gsl_vector_get(g2->subDomain[l2][m2][n2].atom[b]->pos,dim);
										}
										if (p->boundary[0]=='p'){
											if (l2>l+1){
												p2[0]-=p->sideL[0];
											}
											if (l2<l-1){
												p2[0]+=p->sideL[0];
											}
										}
										if (p->boundary[1]=='p'){
											if (m2>m+1){
												p2[1]-=p->sideL[1];
											}
											if (m2<m-1){
												p2[1]+=p->sideL[1];
											}
										}
										if (p->boundary[2]=='p'){
											if (n2>n+1){
												p2[2]-=p->sideL[2];
											}
											if (n2<n-1){
												p2[2]+=p->sideL[2];
											}
										}
										for (int dim=0;dim<3&&inBounds;dim++){
											inBounds=((p2[dim]>=outerBound[dim][0])&&(p2[dim]<=outerBound[dim][1]));
											if (inBounds){
												if (p2[dim]>innerBound[dim][0]&&p2[dim]<innerBound[dim][1]){
													nInner++;
												}
											}
										}

										if (inBounds &&nInner==3){
											inBounds=false;
										}

										if (inBounds){
											// Positions already corrected for PBCs
											double dSq=0;

											gsl_vector_view p_2=gsl_vector_view_array(&p2[0],3);
											gsl_vector_sub(&p_2.vector, subDomain[i][j][k].atom[a]->pos);
											gsl_blas_ddot(&p_2.vector, &p_2.vector, &dSq);

											if (dSq<=RcSq&&dSq>=R0Sq){
												tempNN[tempn_NN]=g2->subDomain[l2][m2][n2].atom[b];
												tempn_NN++;
											}
										}
									}
								}
							}
						}
					}
					//Copy NN list to atoms
					subDomain[i][j][k].atom[a]->n_NN=tempn_NN;
					AllocAtomNN(subDomain[i][j][k].atom[a]);
					for(int l=0;l<tempn_NN;l++){
						subDomain[i][j][k].atom[a]->NN[l]=tempNN[l];
					}
				}
			}
		}
	}

	gsl_vector_free(d);
	delete [] tempNN;
}
void Group::CalcAvgNeighDistVsR(Distrib *NN_r)
{
	//pos mags for atoms in group in GetMaxR(c) command during alloc of NN_r
	double *n=new double[NN_r->n];
	gsl_vector *dist=gsl_vector_calloc(3);
	for(int i=0;i<NN_r->n;i++){
		n[i]=0;
	}
	double L_2[3];
	for (int i=0;i<3;i++){
		L_2[i]=p->sideL[i]/2.0;
	}
	for (int i=0;i<nAtoms;i++){
		double d=0;
		int k=(int)(atom[i]->r/NN_r->step);
		for (int j=0;j<atom[i]->n_NN;j++){
			CalcDistPBC_gsl(atom[i]->NN[j]->pos, atom[i]->pos, dist, &p->boundary[0], &L_2[0], &p->sideL[0]);
			double temp;
			gsl_blas_ddot(dist, dist, &temp);
			d+=sqrt(temp);
		}
		//Keep seperate running totals and divide in the end
		NN_r->y[k]+=d;
		n[k]+=atom[i]->n_NN;
	}
	double dTot=0, nD=0;
	for (int i=0;i<NN_r->n;i++){
		if (n[i]!=0){
			dTot+=NN_r->y[i];
			nD+=n[i];
			NN_r->y[i]/=n[i];
		}
	}
	cout<<"Average Neighbor Distance: "<<dTot/nD<<endl;
	delete [] n;
	gsl_vector_free(dist);

}
void Group::CalcNumNNDist(Distrib *nNN, double maxNN)
{
	int max=(int) maxNN;
	for (int i=0;i<nAtoms;i++){
		if (atom[i]->n_NN<=(int) max){
			nNN->y[atom[i]->n_NN]++;
		}
	}
}
void Group::AllocSubDomains(int nX, int nY, int nZ)
{
	try{
		if (subDomain){
			for(int i=0;i<nSubDomX;i++){
				for (int j=0; j<nSubDomY; j++){
					for (int k=0;k<nSubDomZ;k++){
						delete [] subDomain[i][j][k].atom;
						subDomain[i][j][k].nAtoms=0;
					}
					//TODO Currently some problem when NN Alloc is called twice in the same run
					delete [] subDomain[i][j];
				}
				delete [] subDomain[i];
			}
			delete [] subDomain;
		}

		nSubDomX=nX;
		nSubDomY=nY;
		nSubDomZ=nZ;

		subDomain=new Group ** [nSubDomX];
		for (int i=0;i<nSubDomX;i++){
			subDomain[i]=new Group * [nSubDomY];
			for (int j=0;j<nSubDomY;j++ ){
				subDomain[i][j]=new Group [nSubDomZ];
			}
		}
	}
	catch(exception &e){
		cout<<"Error in SubDomain Alloc: "<<e.what()<<endl;
		exit(1);
	}
}

void Group::CreateSubDomains(int *nSD)
{
	AllocSubDomains(nSD[0], nSD[1],nSD[2]);
	double Ls[3];
	for (int i=0;i<3;i++){
		Ls[i]=p->sideL[i]/(double) nSD[i];
	}
	int **address=new int *[nAtoms];
	for (int i=0;i<nAtoms;i++){
		address[i]=new int [3];
	}

	for (int i=0;i<nAtoms;i++){
		for (int j=0;j<3;j++){
			address[i][j]=(int)((atom[i]->pos->data[j]-p->dim[j][0])/Ls[j]);
		}
		// For some reason, in the release version address can be found to be equal to the dimension of the subdomains.
		// This should not be because the atoms are already checked for this in function Group::RemapPos() ...
		// The following is a further check and adds these atom into the box near the surface.
		if (address[i][0]==nSubDomX){
			address[i][0]--;
		}
		if (address[i][1]==nSubDomY){
			address[i][1]--;
		}
		if (address[i][2]==nSubDomZ){
			address[i][2]--;
		}

		subDomain[address[i][0]][address[i][1]][address[i][2]].nAtoms++;
	}
	for (int i=0;i<nSubDomX;i++){
		for(int j=0;j<nSubDomY;j++){
			for(int k=0;k<nSubDomZ;k++){
				try{
					subDomain[i][j][k].atom=new Atom * [subDomain[i][j][k].nAtoms];
					subDomain[i][j][k].nAtoms=0;
				}
				catch(exception &e){
					cout<<"Exception3: "<<e.what()<<endl;
				}
			}
		}
	}
	int a,b,c;
	for(int i=0;i<nAtoms;i++){
		a=address[i][0];
		b=address[i][1];
		c=address[i][2];
		try{
			subDomain[a][b][c].atom[subDomain[a][b][c].nAtoms]=atom[i];
			subDomain[a][b][c].nAtoms++;
		}
		catch(exception &e){
			cout<<"Exception2: "<<e.what()<<endl;
		}
	}

	for (int i=0;i<nAtoms;i++){
		delete [] address[i];
	}
	delete [] address;
}
void Group::OutputPosition(string format, string path, string name)
{
	//Default output variables
	bool scalePos=false;
	bool writeIndex=true;
	bool writeAtomType=true;
	string fileName=path+"/"+name;
	ofstream ofile;
	ofile.setf(ios::fixed, ios::floatfield);
	//Write header to file depending on file type
	//also define variables for format position output.
	if (format=="lammpstrj"){
		scalePos=true;
		writeIndex=true;
		writeAtomType=true;
		ofile.open(fileName.c_str(), ios_base::app);
		if (ofile.is_open()){
			ofile<<"ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n";
			ofile<<nAtoms<<endl;
			ofile<<"ITEM: BOX BOUNDS\n";
			for (int i=0;i<3;i++){
				ofile<<p->dim[i][0]<<" "<<p->dim[i][1]<<endl;
			}
			ofile<<"ITEM: ATOMS\n";
			ofile.close();
		}

	}
	else{
		cout<<"Output of average position in format: "<<format<<" is not supported.\n";
		exit(0);
	}
	ofile.open(fileName.c_str(), ios_base::app);
	if (ofile.is_open()){
		for (int i=0;i<nAtoms;i++){
			if (writeIndex){
				ofile<<atom[i]->index<<" ";
			}
			if (writeAtomType){
				ofile<<atom[i]->type<<" ";
			}
			for (int j=0;j<3;j++){
				if (scalePos){
					ofile<<(atom[i]->pos->data[j]-p->dim[j][0])/p->sideL[j]<<" ";
				}
				else{
					ofile<<atom[i]->pos->data[j]<<" ";
				}
			}
			ofile<<endl;
		}
		ofile.close();
	}

	cout<<"Output position for group "<<this->name<<" to file: "<<fileName<<endl;
}
bool Group::CheckAtomsInBox()
{
	for (int i=0;i<nAtoms;i++){
		for (int j=0;j<3;j++){
			if (atom[i]->pos->data[j]<p->dim[j][0]||atom[i]->pos->data[j]>p->dim[j][1]){
				return false;
			}
		}
	}
	return true;
}

void Group::ShrinkGroup()
{
	Atom **tempAtoms=new Atom *[nAtoms];
	int nonZeroAtoms=0;
	for (int i=0;i<nAtoms;i++){
		if (atom[i]!=0){
			tempAtoms[nonZeroAtoms]=atom[i];
			nonZeroAtoms++;
		}
	}
	delete [] atom;
	nAtoms=nonZeroAtoms;
	atom=new Atom *[nAtoms];
	for(int i=0; i<nAtoms; i++){
		atom[i]=tempAtoms[i];
	}
	delete [] tempAtoms;
}

void Group::NoGroupBoolean(string booleanName)
{
	cout<<"Error: Cannot create group: "<<name<<", group creation boolean: "<<booleanName<<" not supported.\n";
	exit(0);
}

double Group::CalcRadiusOfGyration(gsl_vector *center)
{
	double rgSq=0;
	gsl_vector *pos=gsl_vector_alloc(3);
	for (int i=0;i<nAtoms;i++){
		double temp=0;
		gsl_vector_memcpy(pos,atom[i]->pos);
		gsl_vector_sub(pos,center);
		gsl_blas_ddot(pos,pos,&temp);
		rgSq+=temp;
	}
	rgSq/=(double)nAtoms;
	gsl_vector_free(pos);
	return sqrt(rgSq);
}
void Group::CalcDensityVsR(Distrib *density)
{
	for (int i=0;i<nAtoms;i++){
		if (atom[i]->r>=density->minX&&atom[i]->r<=density->maxX){
			density->y[(int) (atom[i]->r/density->step)]++;
		}
	}
	for (int i=0;i<density->n;i++){
		if (density->x[i]!=0){
			density->y[i]/=4*PI*density->x[i]*density->x[i]*density->step;
		}
	}
}
