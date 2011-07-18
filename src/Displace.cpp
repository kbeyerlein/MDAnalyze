/*
 * Displace.cpp
 *
 *  Created on: Jan 11, 2010
 *      Author: Ken
 */

#include "Displace.h"

Displace::Displace() {
	pos_i=0;
	pos_f=0;
	nAtoms=0;
	u=0;
	uSq=0;
	specifyMode="";
	groupName="";
	initial=-1;
	final=-1;
	radMSD=0;
	pMSD=0;
	radDisp=0;
	msd=0;
	scaled=false;
}

Displace::~Displace() {
	if (u){
		for (int i=0; i<nAtoms;i++){
			gsl_vector_free( u[i]);
		}
		delete [] u;
	}
	delete [] uSq;
	if(pMSD){
		CleanDistrib(pMSD);
		delete pMSD;
	}
	if(radMSD){
		CleanDistrib(radMSD);
		delete radMSD;
	}
	if(radDisp){
		CleanDistrib(radDisp);
		delete radDisp;
	}
}

void Displace::Allocate()
{
	try{
		if (u){
			for (int i=0;i<nAtoms;i++){
				if (u[i]){
					gsl_vector_free(u[i]);
				}
			}
			delete [] u;
		}
		if (uSq){
			delete [] uSq;
		}

		if(pos_i->nAtoms==pos_f->nAtoms){
			nAtoms=pos_i->nAtoms;
		}
		else{
			cout<<"Cannot calc displacements between time steps"<<endl;
			cout<<"	Not the same number of atoms"<<endl;
		}
		u=new gsl_vector* [nAtoms];
		for (int i=0;i<nAtoms;i++){
			u[i]=gsl_vector_calloc(3);
		}
		uSq=new double [nAtoms];
	}
	catch(exception &e){
		cout<<"Exception in Displace::Allocate(): "<<e.what();
		exit(1);
	}
	Zero();
}
void Displace::CalcDisplacements(double scale_i, double scale_f)
{
	if(u!=0){
		pos_i->CalcCenter();
		pos_f->CalcCenter();
		for (int i=0;i<nAtoms;i++){
			if (pos_i->atom[i]->index==pos_f->atom[i]->index){
				gsl_vector *temp_i=gsl_vector_alloc(3);
				gsl_vector *temp_f=gsl_vector_alloc(3);
				gsl_vector_memcpy(temp_i, pos_i->atom[i]->pos);
				gsl_vector_memcpy(temp_f, pos_f->atom[i]->pos);
				// Use positions from center of mass
				gsl_vector_sub(temp_i, pos_i->p->center);
				gsl_vector_sub(temp_f, pos_f->p->center);
				// Scale vectors as desired
				gsl_vector_scale(temp_i, scale_i);
				gsl_vector_scale(temp_f, scale_f);
				// Calculated displacement vector, u
				gsl_vector_sub(temp_f,temp_i);
				gsl_vector_memcpy(u[i],temp_f);
				// Calc uSq
				gsl_blas_ddot(u[i],u[i], &(uSq[i]));
				gsl_vector_free(temp_i);
				gsl_vector_free(temp_f);
			}
			else{
				cout<<"Cannot calculate displacement of atom: "<<pos_i->atom[i]->index<<endl;
			}
		}
	}
	else{
		cout<<"Error: Displacement array has not been allocated.\n";
		exit(0);
	}
}
double Displace::GetMaxDisplaceSq()
{
	double max=uSq[0];
	for (int i=1;i<nAtoms;i++){
		max =getmax(max,uSq[i]);
	}
	if (max==0){
		cout<<"Warning max uSq=0, Likely that displacements have not been calculated.\n";
	}
	return max;
}
void Displace::CalcPMSD()
{
	AllocDistrib(&pMSD, "pMSD", "Displacement", "p(MSD)", 0, GetMaxDisplaceSq(), DPREC);
	double increment=1/(double)nAtoms;
	for (int i=0;i<nAtoms;i++){
		pMSD->y[(int)(uSq[i]/pMSD->step)]+=increment;
	}
}

void Displace::CalcRadMSD(gsl_vector *center)
{
	//Center is calculated outside of routine, seperately
	pos_i->CreateRadialDistrib(center, DPREC*10);
	pos_i->CalcPosMags(center);
	AllocDistrib(&radMSD, "RadMSD", "RadialPosition", "AvgMSD", 0, pos_i->GetMaxR(), DPREC*10);
	for (int i=0;i<nAtoms;i++){
		radMSD->y[(int)(pos_i->atom[i]->r/radMSD->step)]+=uSq[i];
	}
	for (int i=0;i<radMSD->n;i++){
		if(radMSD->y[i]!=0){
			if (pos_i->p->rad->y[i]!=0){
				radMSD->y[i]/=pos_i->p->rad->y[i];
			}
			else{
				cout<<"Error: radial distribution and radMSD do not match up! "<<i<<"\n";
				exit(0);
			}
		}
	}
}
void Displace::CalcRadDisp(gsl_vector *center)
{
	//Center is calculated outside of routine, seperately
	pos_i->CreateRadialDistrib(center, DPREC*10);
	pos_i->CalcPosMags(center);
	AllocDistrib(&radDisp, "RadDisp", "RadialPosition", "AvgRadDisplacement", 0, pos_i->GetMaxR(), DPREC*10);
	gsl_vector *centPos=gsl_vector_alloc(3);
	for (int i=0;i<nAtoms;i++){
		double del=0;
		// Create centered position vector
		gsl_vector_memcpy(centPos,pos_i->atom[i]->pos);
		gsl_vector_sub(centPos, center);

		gsl_blas_ddot(u[i], centPos, &del);
		if (pos_i->atom[i]->r>=pos_i->p->rad->step){
			del/=pos_i->atom[i]->r;
		}
		else{
			del=0;
		}
		radDisp->y[(int)(pos_i->atom[i]->r/radDisp->step)]+=del;

	}
	gsl_vector_free(centPos);
	for (int i=0;i<radDisp->n;i++){
		if (radDisp->y[i]!=0){
			if (pos_i->p->rad->x[i]==radDisp->x[i]){
				radDisp->y[i]/=(double)pos_i->p->rad->y[i];
			}
			else{
				cout<<"Error: radial distribution and radDisp do not match up! "<<i<<endl;
				exit(0);
			}
		}
	}
}
/*
void Displace::CalcScaledDisplacement(double scale_i, double scale_f)
{
	if(u!=0){
		pos_i->CalcCenterOfMass();
		pos_f->CalcCenterOfMass();
		for (int i=0;i<nAtoms;i++){
			if (pos_i->atom[i].index==pos_f->atom[i].index){
				u_2[i]=0;
				for (int j=0;j<3;j++){
					u[i][j]=(pos_f->atom[i].pos[j]-pos_f->center[j])/scale_f-(pos_i->atom[i].pos[j]-pos_i->center[j])/scale_i;
					u_2[i]+=u[i][j]*u[i][j];
				}
			}
			else{
				cout<<"Cannot calculate displacement of atom: "<<pos_i->atom[i].index<<endl;
			}
		}
	}
	else{
		cout<<"Error: Displacement array has not been allocated.\n";
		exit(0);
	}
}
*/
void Displace::CalcAvgMSD()
{
	msd=0;
	for (int i=0;i<nAtoms;i++){
		msd+=uSq[i];
	}
	msd/=(double)nAtoms;
	if (msd==0){
		cout<<"Warning Avg MSD=0, Likely that displacements have not been calculated.\n";
	}
}
/*
void Displace::OutputMSD(string path, string _name)
{
	string _file=path+"/"+_name+"MSD.lammpstrj";
	double min=u_2[0], max=u_2[0];
	for (int i=0;i<nAtoms;i++){
		min=getmin(min,u_2[i]);
		max=getmax(max,u_2[i]);
	}
	double nInd=10;
	double step=(max-min)/nInd;
	cout<< "\nOutputting MSD/atom to File: "<<_file<<endl;
	fstream file(_file.c_str(), fstream::out|fstream::app);

	bool scale=pos_i->ScalePosition();

	file << "ITEM: TIMESTEP\n"<<pos_i->timeStep<<"\n" << "ITEM: NUMBER OF ATOMS\n" << nAtoms << "\n";
	file << "ITEM: BOX BOUNDS\n" ;
	double L[3];
	for (int i=0;i<3;i++){
		for(int j=0;j<2;j++){
			file<< pos_i->dim[i][j] << " ";
		}
		L[i]=pos_i->dim[i][1]-pos_i->dim[i][0];
		file<<endl;
	}
	file << "ITEM: ATOMS\n";
	for (int i=0;i<pos_i->nAtoms;i++){
		file.precision(0);
		file << pos_i->atom[i].index << " " << (int)((u_2[i]-min)/step)<< " ";
		file.precision(5);
		for (int k=0;k<3;k++){
			if (scale){
				file<< ((pos_i->atom[i].pos[k]-pos_i->center[k])/L[k])+0.5<< " " ;
			}
			else{
				file<< pos_i->atom[i].pos[k]<< " " ;
			}
		}
		file<<endl;
	}
	file.close();
}
*/
void Displace::Zero()
{
	for (int i=0;i<nAtoms;i++){
		uSq[i]=0;
	}
}
