/*
 * Position.cpp
 *
 *  Created on: Jan 3, 2010
 *      Author: Ken
 */
#include "Position.h"
Position::Position()
{
	file="";
	fileType="";
	rad=0;
	rdf=0;
	center=0;
	for (int i=0;i<3;i++){
		for(int j=0;j<2;j++){
			dim[i][j]=0;
		}
		sideL[i]=0;
		boundary[i]='p';
	}
}

Position::~Position(void)
{
	if (center){
		gsl_vector_free(center);
	}
	if(rad){
		CleanDistrib(rad);
		delete rad;
	}
	if(rdf){
		CleanDistrib(rdf);
		delete rdf;
	}
	if (lattice.points){
		for (int i=0;i<lattice.nPoints;i++){
			delete [] lattice.points[i];
		}
		delete [] lattice.points;
	}
}
/*
 //TODO Write Input routines for other position file types
void Position::ReadPosFile()
{
	string tempstring1, tempstring2;
	bool stop=false;
	ifstream infile;
	infile.open(fileName.c_str());
	if (infile.is_open()){
		printf("file %s opened\n", fileName.c_str());
		while ((!infile.eof())&& !stop)
		{
			infile >> tempstring1;

			//Input number of atoms in simulation
			if (tempstring1=="NUMBER"){
				infile >> tempstring2>>tempstring2>>nAtoms;
				//allocate mem for positions and other info
				Allocate();
			}
			//Input atom positions
			else if (tempstring1=="ATOMS"){
				for (int i=0; i<nAtoms; i++){
					infile >>atom[i].index>>atom[i].type>>atom[i].pos[0]>> atom[i].pos[1] >>atom[i].pos[2];
					atom[i].r=-1;
					atom[i].n_NN=0;
					atom[i].NN=0;

				}
				stop=true;
			}
		}
		infile.close();
	}
	else {
		printf("Could not open file: %s \n", fileName.c_str());
		exit (0);
	}
}
*/

/*
void Position::ReadLammpstrjFile(string path, string file)
{
	//string tempstring1, tempstring2;
	fileName=path+"/"+file;
	string tempstring, item;
	bool stop=false, input=false;
	double sideL[3];
	ifstream infile;
	infile.open(fileName.c_str());
	if (infile.is_open()){
		printf("file %s opened\n", fileName.c_str());
		while ((!infile.eof())&& !stop)
		{
			infile >> tempstring >>item;

			//Input number of atoms in simulation
			if (item=="TIMESTEP"){
				int step;
				infile>>step;
				if (timeStep==step){
					input=true;
				}
			}
			else if(item=="NUMBER"){
				infile.ignore(100,'\n');
				infile >> nAtoms;
				//allocate mem for positions and other info
				Allocate();
			}
			else if (item=="BOX"){
				if (input){
					infile.ignore(100,'\n');
					for (int i=0;i<3;i++){
						for(int j=0;j<2;j++){
							infile>>dim[i][j];
						}
						sideL[i]=dim[i][1]-dim[i][0];
					}
				}
				else{
					for (int i=0;i<4;i++){
						infile.ignore(100,'\n');
					}
				}
			}
			//Input atom positions
			else if (item=="ATOMS"){
				if (input){
					infile.ignore(100,'\n');
					for (int i=0; i<nAtoms; i++){
						infile >>atom[i].index>>atom[i].type>>atom[i].pos[0]>> atom[i].pos[1] >>atom[i].pos[2];

					}
					bool scale=ScalePosition();
					for (int i=0;i<nAtoms;i++){
						atom[i].r=-1;
						atom[i].n_NN=0;
						atom[i].NN=0;
						if (scale){
							for (int j=0;j<3;j++){
								atom[i].pos[j]*=sideL[j];
								atom[i].pos[j]+=dim[j][0];
							}
						}
					}
					stop =true;
				}
				else{
					for (int i=0;i<nAtoms+1;i++){
						infile.ignore(250,'\n');
					}
				}
			}
			else{
				cout<<"Unsupported keyword in lammpstrj file: "<<item<<endl;
				exit(0);
			}
		}
		if (!stop){
			cout<<"Could not find timestep: "<<timeStep<<" in lammpstrj file."<<endl;
			exit(0);
		}
		infile.close();
	}
	else {
		printf("Could not open file: %s \n", fileName.c_str());
		exit(0);
	}
}


int Position::ReadMoPosFile()
{
	int count=0;
	ifstream infile;
	infile.open(fileName.c_str());
	if (infile.is_open()){
		printf("file %s opened\n", fileName.c_str());
		while ((!infile.eof()))
		{
			infile.ignore(256, '\n');
			count++;
		}
		infile.close();
	}
	else {
		printf("Could not open file: %s \n", fileName.c_str());
		return -1;
	}
	nAtoms=count;
	x=(float*) _mm_malloc(nAtoms*sizeof(float), 16);
	y=(float*) _mm_malloc(nAtoms*sizeof(float), 16);
	z=(float*) _mm_malloc(nAtoms*sizeof(float), 16);
	if (x==0||y==0||z)
	infile.open(fileName.c_str());
	if (infile.is_open()){
		for (int i=0; i<nAtoms; i++){
			infile >>x[i]>> y[i] >>z[i];
			infile.ignore(256, '\n');

		}
		infile.close();
	}
	return 0;
}
int Position::ReadTaoPosFile()
{
	int count=0;
	ifstream infile;
	infile.open(fileName.c_str());
	if (infile.is_open()){
		printf("file %s opened\n", fileName.c_str());
		infile>>nAtoms;
		for (int i=0;i<3;i++){
			for (int j=0;j<2;j++){
				infile>>dim[i][j];
			}
		}
		x=(float*) _mm_malloc(nAtoms*sizeof(float), 16);
		y=(float*) _mm_malloc(nAtoms*sizeof(float), 16);
		z=(float*) _mm_malloc(nAtoms*sizeof(float), 16);
		int i;
		for (i=0; i<nAtoms&&!infile.eof(); i++){
			infile >>x[i]>> y[i] >>z[i];
		}
		float temp;
		while(!infile.eof()){
			infile>>temp;
			count++;
		}
		if (nAtoms!=i||count>1){
			printf("Inaccurate number of atoms specified in Tao's position file. Double Check!!\n");
			return -1;
		}
		infile.close();
	}
	GetDimensions();
	return 0;
}
*/

double Position::GetMaxD()
{
	double maxX=0;
	for (int i=0;i<3; i++){
		double temp=dim[i][1]-dim[i][0];
		maxX+=temp*temp;
	}
	return sqrt(maxX);
}
//TODO Rewrite functions for new structures

/*
void Position::CalcLatticeVectors()
{
	double **perfectNNvect, *perfectNNDistSq;
	Lattice *localStruct;
	int nNNPerfVect;
	if (lattice.name=="primative"){
		nNNPerfVect=3;
	}
	else{
		nNNPerfVect=lattice.nPoints-1;
	}
	try{
		perfectNNvect=new double *[nNNPerfVect];
		for (int i=0;i<nNNPerfVect;i++){
			perfectNNvect[i]=new double [3];
		}
		perfectNNDistSq=new double [nNNPerfVect];
		localStruct=new Lattice [nAtoms];
	}
	catch(exception &e){
		cout<<"Error in CalcLatticeVectors: " <<e.what();
		exit(1);
	}
	//Create perfect primary lattice NN
	for(int i=0;i<nNNPerfVect;i++){
		for (int dim=0;dim<3;dim++){
			if (lattice.name=="primative"){
				perfectNNvect[i][dim]=lattice.vectors[i][dim];
			}
			else{
				perfectNNvect[i][dim]=0;
				for (int j=0;j<3;j++){
					perfectNNvect[i][dim]+=lattice.points[i+1][j]*lattice.vectors[j][dim];
				}
			}
		}
	}
	for (int i=0;i<nNNPerfVect;i++){
		perfectNNDistSq[i]=Dot(perfectNNvect[i], perfectNNvect[i]);
	}

	int nNN=GetnNN(&lattice);
	for (int i=0;i<nAtoms;i++){
		if (atom[i].n_NN==nNN){
			int *indexNNVect;
			double **currentNNVect;
			try{
				currentNNVect=new double *[nNN];
				for (int j=0;j<nNN;j++){
					currentNNVect[j]=new double [3];
				}
				indexNNVect=new int [nNNPerfVect];
			}
			catch(exception &e){
				cout<<"Error in CalculateLatticeVectors: "<<e.what();
				exit(1);
			}
			//Calc NN vectors for current atom
			for (int j=0;j<nNN;j++){
				for (int dim=0;dim<3;dim++){
					currentNNVect[j][dim]=atom[i].NN[j]->pos[dim]-atom[i].pos[dim];
				}
			}
			//Identify primary NN
			for (int j=0;j<nNNPerfVect;j++){
				indexNNVect[j]=-1;
				double maxDot=0;
				for (int k=0;k<nNN;k++){
					double rSq=Dot(currentNNVect[k],currentNNVect[k]);
					//Check if length is correct for each NN found (more relevant for noncubic lattice)
					if(rSq<1.21*perfectNNDistSq[j]&&rSq>.81*perfectNNDistSq[j]){
						double temp=maxDot;
						//Find NN which is in the closest orientation with the perfect NN
						maxDot=getmax(maxDot,Dot(currentNNVect[k], perfectNNvect[j]));
						if (maxDot>temp){
							indexNNVect[j]=k;
						}
					}
				}
			}
			//Check for errors in identifying primary NN
			bool err=false;
			for (int j=0;j<nNNPerfVect;j++){
				if (indexNNVect[j]==-1){
					cout<<"Error: Could not identify nearest neighbor. Atom: "<<i<<"NN: "<<j<<endl;
					err=true;
				}
				for (int k=j+1;k<nNNPerfVect; k++){
					if (indexNNVect[j]==indexNNVect[k]){
						cout<<"Error: Identified the same NN: "<<j<<" for atom: "<<i<<endl;
						err=true;
					}
				}
			}
			//Can do plane calc of alberto
			//Mine is simpler (maybe not better) where I only consider the primary NN
			if (!err){
				if (lattice.name=="faceCentered"&&lattice.symmetry=="cubic"){

				}
			}

			//Clean up
			for (int i=0;i<nNN; i++){
				delete [] currentNNVect[i];
			}
			delete [] currentNNVect;
			delete [] indexNNVect;
		}
	}
	for (int i=0;i<nNN; i++){
		delete [] perfectNNvect[i];
	}
	delete [] perfectNNvect;
	delete [] perfectNNDistSq;

}
*/
/*
void Position::CalcLatticeVectors()
{
	double refLattNNPlaneL[3];
	lattice.CalcRecipLattVects();
	// Get the length of expected cross product vector to identify correct planes
	lattice.GetCrossNNVectL(&refLattNNPlaneL[0]);
	//Only consider atoms with correct NN for that lattice (ie: FCC=12)
	int nNN=lattice.GetnNN();
	for (int i=0;i<nAtoms;i++){
		if (atom[i].n_NN==nNN){
			//Calc NN Vectors
			gsl_vector **NNVect;
			int *planeIndex;
			try{
				NNVect=(gsl_vector**) malloc(nNN*sizeof(gsl_vector*));
				for(int j=0;j<nNN;j++){
					NNVect[j]=gsl_vector_calloc(3);
				}
				planeIndex=new int [nNN];
			}
			catch(exception &e){
				cout<<"Exception found in CalcLatticeVectors: "<<e.what()<<endl;
				exit(1);
			}
			for (int j=0;j<nNN;j++){
				for (int k=0;k<3;k++){
					gsl_vector_set(NNVect[j],k,atom[i].NN[j]->pos[k]-atom[i].pos[k]);
				}
			}
			//Find Correct Plane (compare with expected cross product above)
			Lattice localLatt;
			FindApproxLattPlanes(NNVect, nNN, refLattNNPlaneL, &localLatt);
			//Determine plane index of each nearest neighbor
			for (int j=0;j<nNN;j++){
				double minDot=0;
				for(int k=0;k<3;k++){
					gsl_vector_view plane=gsl_vector_view_array(localLatt.recipVectors[j],3);
					double dot=0;
					gsl_blas_ddot(NNVect[j],&plane.vector, &dot);
					if (k==0){
						minDot=abs(dot);
						planeIndex[j]=0;
					}
					else{
						if (minDot>abs(dot)){
							minDot=abs(dot);
							planeIndex[j]=k;
						}
					}
				}
			}
			//Determine nAtoms in each plane
			int nAtomsInPlane[3];
			for (int j=0;j<3;j++){
				nAtomsInPlane[j]=1; //Init to one to count atom in center which is on all planes
			}
			for (int j=0;j<nNN;j++){
				nAtomsInPlane[planeIndex[j]]++;
			}
			//Insert routine to find best plane between all atoms (Lagrangian Mult or Avg of all possible)
			double distToPlane[3];
			for (int j=0;j<3;j++){
				gsl_vector **pointsInPlane;
				double *weight;
				try{
					pointsInPlane=(gsl_vector **) malloc(nAtomsInPlane[j]*sizeof(gsl_vector*));
					for (int k=0;k<nAtomsInPlane[j]; k++){
						pointsInPlane[k]=gsl_vector_calloc(3);
					}
					weight=new double [nAtomsInPlane[j]];
				}
				catch(exception &e){
					cout<<"Exception found in CalcLatticeVectors: " <<e.what()<<endl;
					exit(1);
				}
				//Init center atom in plane
				for (int k=0;k<3;k++){
					gsl_vector_set(pointsInPlane[0],k,0);
				}
				int count=1;
				//Copy atoms in the plane
				for (int k=0;k<nNN;k++){
					if (planeIndex[k]==j){
						gsl_vector_memcpy(pointsInPlane[count], NNVect[k]);
						//Assuming identical atoms
						weight[count]=1.0;
						count++;
					}
				}
				//Find the plane which best approximates the set of points
				//Output plane is in terms of the unit normal vector.
				gsl_vector_view plane=gsl_vector_view_array(localLatt.recipVectors[j], 3);
				FindBestPlaneOfPoints(pointsInPlane,weight, nAtomsInPlane[j], &plane.vector, &distToPlane[j]);

				for(int k=0;k<nAtomsInPlane[j];k++){
					gsl_vector_free(pointsInPlane[k]);
				}
				free(pointsInPlane);
				delete [] weight;
			}
			//Generate local real lattice vectors from reciprocal vectors
			//With normalizations loose information on extent of lattice deformation.
			for (int j=0;j<3;j++){
				Normalize(localLatt.recipVectors[j]);
			}
			localLatt.CalcRealLattVects();
			for (int j=0;j<3;j++){
				Normalize(localLatt.vectors[j]);
			}

			//


			for (int j=0;j<nNN;j++){
				gsl_vector_free(NNVect[j]);
			}
			free(NNVect);
			delete [] planeIndex;
		}
	}
}
*/
void Position::FindApproxLattPlanes(gsl_vector **NN, int nNN, double *refLattNNPlaneL, Lattice *localLatt)
{
	bool planeFound[3];
	for (int i=0;i<3;i++){
		planeFound[i]=false;
	}
	bool foundAll=false;
	for (int i=0;i<nNN&&!foundAll;i++){
		gsl_vector* tempB=gsl_vector_calloc(3);
		bool stop=false;
		for (int j=i+1;j<nNN&&!stop;j++){
			Cross(NN[i],NN[j],tempB);
			double l=VectorL(tempB);
			int index;
			for(int m=0;m<3;m++){
				//Allowing 10% variance (arbitrary and should test)
				if (l>.9*refLattNNPlaneL[m]&&l<1.1*refLattNNPlaneL[m]&&!planeFound[m]){
					index=m;
					double maxDot=0;
					gsl_vector* negB=gsl_vector_alloc(3);
					gsl_vector_memcpy(negB, tempB);
					gsl_vector_scale(negB, -1);
					//Check other recip vectors with same length to see if a better match is possible.
					bool useNegVect=false;
					for (int n=m; n<m+3; n++){
						if (refLattNNPlaneL[m]==refLattNNPlaneL[n%3]&&!planeFound[n%3]){
							double temp=maxDot;
							double dot=0;
							gsl_vector_view recip= gsl_vector_view_array(lattice.recipVectors[n%3],3);
							gsl_blas_ddot(tempB,&recip.vector,&dot);
							maxDot=getmax(maxDot, dot);
							gsl_blas_ddot(negB,&recip.vector, &dot);
							maxDot=getmax(maxDot, dot);
							if (maxDot>temp){
								index=n%3;
								if (maxDot==dot){
									useNegVect=true;
								}
							}
						}
					}
					planeFound[index]=true;
					stop=true;
					if (planeFound[(index+1)%3]&&planeFound[(index+2)%3]) foundAll=true;
					//Save Found plane, make sure correct direction is considered
					if (!useNegVect){
						for (int n=0;n<3;n++){
							localLatt->recipVectors[index][n]=gsl_vector_get(tempB,n);
						}
					}
					else{
						for (int n=0;n<3;n++){
							localLatt->recipVectors[index][n]=gsl_vector_get(negB, n);
						}
					}
					gsl_vector_free(negB);
				}

			}
		}
		gsl_vector_free(tempB);
	}
}
bool Position::SameBoxInfo(Position *pos2)
{
	for (int i=0;i<3;i++){
		if (boundary[i]!=pos2->boundary[i]){
			return false;
		}
		if (sideL[i]!=pos2->sideL[i]){
			return false;
		}
		if (dim[i][0]!=pos2->dim[i][0]){
			return false;
		}
		if (dim[i][1]!=pos2->dim[i][1]){
			return false;
		}
	}
	return true;
}

