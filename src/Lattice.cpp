/*
 * Lattice.cpp
 *
 *  Created on: Mar 17, 2010
 *      Author: Ken
 */

#include "Lattice.h"

Lattice::Lattice() {
	name="";
	symmetry="";
	a=0;
	b=0;
	c=0;
	alpha=0;
	beta=0;
	gamma=0;
	for (int i=0;i<3;i++){
		for (int j=0;j<3;j++){
			vectors[i][j]=0;
			recipVectors[i][j]=0;
		}
	}
	nPoints=0;
	points=0;

}

Lattice::~Lattice() {
	if (points){
		for (int i=0;i<nPoints;i++){
			delete [] points[i];
		}
		delete [] points;
	}
}
double Lattice::GetNNDist()
{
	if (name=="faceCentered"&&symmetry=="cubic"){
		return a/sqrt(2.0);
	}
	else {
		cout<<"Lattice: "<<name<<" and symmetry: "<<symmetry<<" not supported in GetNNDist.\n";
		exit(0);
	}
}
double Lattice::GetNN2Dist()
{
	if (name=="faceCentered"&&symmetry=="cubic"){
		return a;
	}
	else{
		cout<<"Lattice: "<<name<<" and symmetry: "<<symmetry<<" not supported in GetNN2Dist.\n";
		exit(0);
	}
}
int Lattice::GetnNN()
{
	if (name=="faceCentered"&&symmetry=="cubic"){
		return 12;
	}
	else{
		cout<<"Lattice: "<<name<<" and symmetry: "<<symmetry<<" not supported in GetnNN.\n";
		exit(0);
	}
}
void Lattice::LoadLatticePoints(string path, string fileName)
{
	bool latticeLoaded=false;
	string temp;
	string file=path+"/"+fileName;
	printf("Reading Lattice Information file: %s\n",file.c_str());
	ifstream input;
	input.open(file.c_str(), ios::in);
	if(input.is_open()){
		while(!input.eof()&&!latticeLoaded){
			input >> temp;
			if (temp==name){

				input >> nPoints;
				//Init latt points (ie: faceCentered has 4 lattice points)
				AllocPoints();
				//Input Points
				for (int i=0;i<nPoints;i++){
					input>>temp;
					if(temp=="{"){
						for (int j=0;j<3;j++){
							input>>points[i][j];
						}
						input>>temp;
						if (temp!="}"){
							cout<<"Lattice "<<name<<" is not specifed correctly... check latticeFile.\n";
							exit(0);
						}
					}
					else{
						cout<<"Lattice "<<name<<" is not specifed correctly... check latticeFile.\n";
						exit(0);
					}
				}
				latticeLoaded=true;
			}
			else {
				input.ignore(256, '\n');
			}
		}
		input.close();
	}
	else{
		cout<<"Could not open Lattice Information File....check the pathname."<<endl;
		exit(0);
	}
	if (!latticeLoaded){
		cout<<"Lattice: "<<name<< " is not supported in lattice file."<<endl;
		exit(0);
	}
}
void Lattice::InitLatticeVect()
{
	vectors[0][0]=a*sin(alpha*PI/180.0);
	vectors[0][1]=0;
	vectors[0][2]=a*cos(alpha*PI/180.0);
	vectors[1][0]=b*sin(beta*PI/180.0)*cos(gamma*PI/180.0);
	vectors[1][1]=b*sin(beta*PI/180.0)*sin(gamma*PI/180.0);
	vectors[1][2]=b*cos(beta*PI/180.0);
	vectors[2][0]=0;
	vectors[2][1]=0;
	vectors[2][2]=c;
}
void Lattice::CalcRecipLattVects()
{
	for (int i=0;i<3;i++){
		Cross(vectors[(i+1)%3],vectors[(i+2)%3],recipVectors[i]);
		ScaleVect(recipVectors[i],1.0/CalcVol("real"));
	}

}
void Lattice::CalcRealLattVects()
{
	for (int i=0;i<3;i++){
		Cross(recipVectors[(i+1)%3],recipVectors[(i+2)%3],vectors[i]);
		ScaleVect(vectors[i],1.0/CalcVol("reciprocal"));
	}
}

void Lattice::AllocPoints()
{
	try{
		if (points){
			for (int i=0;i< nPoints;i++){
				delete [] points[i];
			}
			delete [] points;
		}
		points= new double *[nPoints];
		for (int i=0;i<nPoints;i++){
			points[i]=new double [3];
		}
	}
	catch(exception &e){
		cout<<"Exception in Lattice::AllocPoints: "<<e.what()<<endl;
		exit(1);
	}
}
void Lattice::CopyLattice(Lattice *to)
{
	to->a=a;
	to->b=b;
	to->c=c;
	to->alpha=alpha;
	to->beta=beta;
	to->gamma=gamma;
	to->name=name;
	to->symmetry=symmetry;
	to->nPoints=nPoints;
	to->AllocPoints();
	for (int i=0;i<to->nPoints;i++){
		for (int j=0;j<3;j++){
			to->points[i][j]=points[i][j];
		}
	}
	for (int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			to->vectors[i][j]=vectors[i][j];
			to->recipVectors[i][j]=recipVectors[i][j];
		}
	}
}
void Lattice::GetCrossNNVectL(double *out)
{
	if (name=="faceCentered"&&symmetry=="cubic"){
		for (int i=0;i<3;i++){
			out[i]=a/2.0;
		}
	}
	else{
		cout<<"Lattice name: " <<name<<" and symmetry: "<<symmetry<<" are not currently supported in Lattice::GetNNRecipVectL().\n";
		exit(0);
	}
}
double Lattice::CalcVol(string space)
{
	if (space=="real"){
		double temp[3];
		Cross(vectors[1], vectors[2], &temp[0]);
		return Dot(vectors[0],&temp[0]);
	}
	else if (space=="reciprocal"){
		double temp[3];
		Cross(recipVectors[1], recipVectors[2], &temp[0]);
		return Dot(recipVectors[0],&temp[0]);
	}
	else{
		cout<<"Cannot calculate volume for space: " << space<<" must be 'real' or 'reciprocal'.\n";
		exit(0);
	}
}
