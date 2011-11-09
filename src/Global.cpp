/*
 * Global.cpp
 *
 *  Created on: Mar 15, 2010
 *      Author: Ken
 */

#include "Global.h"


double Dot(double *x1, double *x2)
{
	return (x1[0]*x2[0]+x1[1]*x2[1]+x1[2]*x2[2]);
}
void Cross(double *x1, double *x2, double *y)
{
	for (int i=0;i<3;i++){
		y[i]=x1[(i+1)%3]*x2[(i+2)%3]-x1[(i+2)%3]*x2[(i+1)%3];
	}
}
void Cross(gsl_vector *x, gsl_vector *y, gsl_vector *z)
{
	for (int i=0;i<3;i++){
		double temp=gsl_vector_get(x,(i+1)%3)*gsl_vector_get(y,(i+2)%3)-gsl_vector_get(x,(i+2)%3)*gsl_vector_get(y, (i+1)%3);
		gsl_vector_set(z,i,temp);
	}
}
void InitBasis(Basis *basis)
{
	bool outUnitCell=false;
	for (int i=0;i<basis->nAtoms;i++){
		for (int j=0;j<3;j++){
			if (gsl_vector_get(basis->atoms[i].pos,i)>1){
				outUnitCell=true;
			}
		}
	}
	if (outUnitCell){
		cout<<"Warning: specified basis places atoms outside of unitcell.\n";
	}
	return;
}
void Normalize(double *x)
{
	double l=sqrt(Dot(x,x));
	if (l!=0){
		ScaleVect(x, 1.0/l);
	}
	else{
		cout<<"Cannot normalize vector, length=0\n";
		exit(0);
	}
}
void ScaleVect(double *x, double scale)
{
	try{
		for (int i=0;i<3;i++){
			x[i]*=scale;
		}
	}
	catch(exception &e){
		cout<<"Exception found in ScaleVect: "<<e.what()<<endl;
		exit(0);
	}
}
double VectorL(gsl_vector *x)
{
	double l;
	gsl_blas_ddot(x,x,&l);
	return sqrt(l);
}
void CalcCenterOfMass_gsl(gsl_vector **points, int n, double *weight, gsl_vector *center)
{
	for (int i=0;i<3;i++){
		gsl_vector_set(center,i,0);
	}
	gsl_vector *wPoints=gsl_vector_alloc(3);
	double avgWeight=0;
	for (int i=0;i<n;i++){
		gsl_vector_memcpy(wPoints,points[i]);
		gsl_vector_scale(wPoints,weight[i]/(double)n);
		gsl_vector_add(center, wPoints);
		avgWeight+=weight[i]/(double)n;
	}
	gsl_vector_scale(center, 1.0/avgWeight);
}
void OuterProduct_gsl(gsl_vector *x, gsl_vector *y, gsl_matrix *result)
{
	for (int i=0;i<(int)x->size;i++){
		double *a=gsl_vector_ptr(x,i);
		for(int j=0; j<(int)y->size; j++){
			double *b=gsl_vector_ptr(y,j);
			gsl_matrix_set(result, i,j,(*a)*(*b));
		}
	}
}
double CalcJacobian_gsl(gsl_matrix *in)
{
	double J;
	gsl_matrix *x=gsl_matrix_alloc(in->size1, in->size2);
	gsl_matrix_memcpy(x, in);
	gsl_permutation *p=gsl_permutation_alloc(3);
	int s;
	gsl_linalg_LU_decomp(x,p,&s);
	J=gsl_linalg_LU_det(x, s);
	gsl_matrix_free(x);
	gsl_permutation_free(p);
	return J;
}


// global distance calc considering PBCs
// Always takes nearest distance of the three projections in a dimension
void CalcDistPBC_gsl(gsl_vector *p1, gsl_vector *p2, gsl_vector *d, char *b, double *L_2, double *L)
{
	// d=p2-p1
	gsl_vector_memcpy(d,p2);
	gsl_vector_sub(d,p1);
	for (int i=0;i<3;i++){
		//Check periodic boundaries
		if (b[i]=='p'){
			if (d->data[i]<0){
				if (d->data[i]<-L_2[i]){
					d->data[i]+=L[i];
				}
			}
			else{
				if (d->data[i]>L_2[i]){
					d->data[i]-=L[i];
				}
			}
		}
	}
}

void QuickSortImproved(double *x, double *y, int lb, int ub)
{
	 while (lb < ub) {
        int m;

        /* quickly sort short lists */
        if (ub - lb <= 50) {
            SortByInsertion(x, y, lb, ub);
            return;
        }

        m = Partition(x, y, lb, ub);

        /* eliminate tail recursion and */
        /* sort the smallest partition first */
        /* to minimize stack requirements    */
        if (m - lb <= ub - m) {
            QuickSortImproved(x, y,lb, m);
            lb = m + 1;
        } else {
            QuickSortImproved(x, y,m + 1, ub);
            ub = m;
        }
    }
}
void SortByInsertion(double *x, double *y, int lb, int ub) {
    int i, j;

    for (i = lb + 1; i <= ub; i++) {
        double t = x[i];
		double u = y[i];

        /* shift down until insertion point found */
		for (j = i-1; j >= lb && (x[j]> t); j--){
            x[j+1] = x[j];
			y[j+1] = y[j];
		}

        /* insert */
        x[j+1] = t;
		y[j+1] = u;
    }
}

int Partition(double *x, double *y, int lb, int ub) {

    /* select a pivot */
    double pivot = x[(lb+ub)/2];

    /* work from both ends, swapping to keep   */
    /* values less than pivot to the left, and */
    /* values greater than pivot to the right  */
    int i = lb - 1;
    int j = ub + 1;
    while (1) {
        double t;
		double u;

        while ((x[--j]> pivot));
        while ((x[++i]< pivot));
        if (i >= j) break;

        /* swap x[i], x[j] */
        t = x[i];
		u = y[i];
        x[i] = x[j];
		y[i]=y[j];
        x[j] = t;
		y[j]=u;
    }

    return j;
}
void CopyAtom(Atom *from, Atom *to)
{
	to->index=from->index;
	to->type=from->type;
	to->r=from->r;
	for(int i=0;i<3;i++){
		if (to->pos&&from->pos){
			gsl_vector_memcpy(to->pos,from->pos);
		}
		if (to->vel&&from->vel){
			gsl_vector_memcpy(to->vel,from->vel);
		}
	}
	to->n_NN=0;
	to->NN=0;
}
void ZeroAtom (Atom *a){
	a->pos=0;
	a->vel=0;
	a->n_NN=0;
	a->index=0;
	a->NN=0;
	a->type=-1;
	a->r=-1;
	a->v=0;
}
void AllocateArray(double **x, int n)
{
	try{
		*x=new double [n];
	}
	catch(exception &e){
		cout<<"Exception in Position::AllocateArray: "<<e.what()<<endl;
		exit(1);
	}
}
void InitDistrib(Distrib *distrib)
{
	for (int i=0;i<distrib->n;i++){
		distrib->x[i]=distrib->minX+i*distrib->step;
		distrib->y[i]=0;
	}
}

void AllocAtomNN(Atom *a)
{
	try{
		if (a->NN) delete [] a->NN;
		a->NN=new Atom *[a->n_NN];
	}
	catch(exception &e){
		cout<<"Error in NN alloc\n";
		exit(1);
	}
}
gsl_vector * ArrayColToGSLVector(double **a, int nX, int y)
{
	gsl_vector *temp=gsl_vector_alloc(nX);
	for (int i=0;i<nX;i++){
		gsl_vector_set(temp,i,a[i][y]);
	}
	return temp;
}
void OutputDistrib(Distrib *dist, string path, string name, bool normX, bool normY, string comment)
{
	double scaleX=dist->x[0], scaleY=dist->y[0];
	if (normX){
		for (int i=1;i<dist->n;i++){
			scaleX=getmax(scaleX,dist->x[i]);
		}
	}
	if (normY){
		for (int i=1;i<dist->n;i++){
			scaleY=getmax(scaleY,dist->y[i]);
		}
	}
	string fileName;
	fileName=path+"/"+name;
	ofstream ofile;
	ofile.setf(ios::scientific, ios::floatfield);
	ofile.open(fileName.c_str(), ios_base::app);
	if (ofile.is_open()){
		if (comment!=""){
			ofile<<comment<<endl;
		}
		ofile.precision(5);
		ofile<<dist->xName<<" ";
		if(normX){
			ofile<<dist->xName+"Norm"<<" ";
		}
		ofile<<dist->yName<<" ";
		if (normY){
			ofile<<dist->yName+"Norm"<<" ";
		}
		ofile<<endl;
		for (int i=0;i<dist->n;i++){
			if(dist->y[i]!=0){
				ofile <<dist->x[i]<<" ";
				if(normX){
					ofile<<dist->x[i]/scaleX<<" ";
				}
				ofile <<dist->y[i]<<" ";
				if (normY){
					ofile<<dist->y[i]/scaleY<<" ";
				}
				ofile<<endl;
			}
		}
		ofile.close();
	}
	else{
		printf("\nCannot open file: %s\n",fileName.c_str());
		exit(0);
	}
	cout<<"Output "<<dist->name<<" file: "<<fileName<<endl;
}
void OutputMultDistribs(Distrib **dist, int nDistribs, string path, string name, bool normX, bool normY, bool sameX, string comment)
{
	double *scaleX=new double [nDistribs], *scaleY=new double [nDistribs];
	if (normX){
		for (int i=0;i<nDistribs;i++){
			if (!sameX||i==0){
				scaleX[i]=dist[i]->x[0];
				for(int j=0;j<dist[i]->n;j++){
					scaleX[i]=getmax(scaleX[i],dist[i]->x[j]);
				}
			}
			else{
				scaleX[i]=scaleX[0];
			}
		}
	}
	if (normY){
		for (int i=0;i<nDistribs;i++){
			scaleY[i]=dist[i]->y[0];
			for(int j=0;j<dist[i]->n;j++){
				scaleY[i]=getmax(scaleY[i],dist[i]->y[j]);
			}
		}
	}
	string fileName;
	fileName=path+"/"+name;
	ofstream ofile;
	ofile.setf(ios::scientific, ios::floatfield);
	ofile.open(fileName.c_str(), ios_base::app);
	if (ofile.is_open()){
		if (comment!=""){
			ofile<<comment<<endl;
		}
		ofile.precision(5);
		for (int i=0;i<nDistribs;i++){
			if (!sameX||i==0){
				ofile<<dist[i]->xName<<" ";
				if(normX){
					ofile<<dist[i]->xName+"Norm"<<" ";
				}
			}
			ofile<<dist[i]->yName<<" ";
			if (normY){
				ofile<<dist[i]->yName+"Norm"<<" ";
			}
		}
		ofile<<endl;
		for (int i=0;i<dist[0]->n;i++){
			bool allZero=true;
			for (int j=0;j<nDistribs;j++){
				if (dist[j]->y[i]!=0){
					allZero=false;
				}
			}
			if (!allZero){
				for (int j=0;j<nDistribs;j++){
					if (!sameX||j==0){
						ofile <<dist[j]->x[i]<<" ";
						if(normX){
							ofile<<dist[j]->x[i]/scaleX[j]<<" ";
						}
					}
					ofile <<dist[j]->y[i]<<" ";
					if (normY){
						ofile<<dist[j]->y[i]/scaleY[j]<<" ";
					}
				}
				ofile<<endl;
			}
		}
		ofile.close();
	}
	else{
		printf("\nCannot open file: %s\n",fileName.c_str());
		exit(0);
	}
	cout<<"Output distributions to file: "<<fileName<<endl;
	delete [] scaleX;
	delete [] scaleY;
}
void AllocDistrib(Distrib **d, string _name, string _xName, string _yName, double _minX, double _maxX, double _step)
{
	try{
		if (*d){
			CleanDistrib(*d);
			delete *d;
		}
		*d=new Distrib;
	}
	catch(exception &e){
		cout<<"Exception found in AllocDistrib: "<<e.what()<<endl;
		exit(1);
	}
	(*d)->name=_name;
	(*d)->xName=_xName;
	(*d)->yName=_yName;
	(*d)->step=_step;
	(*d)->minX=_minX;
	(*d)->maxX=_maxX;
	(*d)->n=(int) (((*d)->maxX-(*d)->minX)/(*d)->step)+1;
	AllocateArray(&(*d)->x, (*d)->n);
	AllocateArray(&(*d)->y, (*d)->n);
	InitDistrib(*d);
}
void AllocDistrib(Distrib *d, string _name, string _xName, string _yName, double _minX, double _maxX, double _step)
{
	(d)->name=_name;
	(d)->xName=_xName;
	(d)->yName=_yName;
	(d)->step=_step;
	(d)->minX=_minX;
	(d)->maxX=_maxX;
	(d)->n=(int) (((d)->maxX-(d)->minX)/(d)->step)+1;
	AllocateArray(&(d)->x, (d)->n);
	AllocateArray(&(d)->y, (d)->n);
	InitDistrib(d);
}
void ScaleDistrib(Distrib *d, string axisName, double scale)
{
	if (axisName==d->xName){
		for (int i=0;i<d->n;i++){
			d->x[i]*=scale;
			d->step*=scale;
		}
	}
	else if (axisName==d->yName){
		for (int i=0;i<d->n;i++){
			d->y[i]*=scale;
		}
	}
	else{
		cout<<"Error: Cannot scale axis: "<<axisName<<" in distribution: "<<d->name<<endl;
		cout<<"       Axis does not exist."<<endl;
		exit(0);
	}
}
void CleanDistrib(Distrib *d)
{
	delete [] d->x;
	delete [] d->y;
}
// Given a set of points this routine finds the best plane through them.
void FindBestPlaneOfPoints(gsl_vector **points, double *weight, int nPoints, gsl_vector *bestPlane, double *distToPlane)
{
	// Calc Center of mass
	gsl_vector *centOfMass=gsl_vector_calloc(3);
	CalcCenterOfMass_gsl(points, nPoints, weight, centOfMass);
	//Translate all points relative to new origin (center of mass)
	for(int i=0;i<nPoints;i++){
		gsl_vector_sub(points[i],centOfMass);
	}
	gsl_matrix *A=gsl_matrix_calloc(3,3);
	//Calc A matrix according to Schomaker (1959) Acta Cryst., 12, 600
	for (int i=0;i<nPoints;i++){
		gsl_matrix *tempA=gsl_matrix_calloc(3,3);
		for (int j=0;j<3;j++){
			double x_j=gsl_vector_get(points[i],j);
			for(int k=0;k<=j;k++){
				double x_k=gsl_vector_get(points[i],k);
				double temp=weight[i]*x_j*x_k/(double)nPoints;
				gsl_matrix_set(tempA,j,k,temp);
				//A is a symmetric matrix
				if (j!=k){
					gsl_matrix_set(tempA,k,j,temp);
				}
			}
		}
		gsl_matrix_add(A,tempA);
		gsl_matrix_free(tempA);
	}
	//alloc variables for eigen calc
	gsl_eigen_symmv_workspace *eigenA=gsl_eigen_symmv_alloc(3);
	gsl_matrix *eigenVectsOfA=gsl_matrix_alloc(3,3);
	gsl_vector *eigenValOfA=gsl_vector_alloc(3);
	//calc eigen vectors and eigen values
	gsl_eigen_symmv(A, eigenValOfA, eigenVectsOfA, eigenA);
	//Eigenvector with smallest eigen value is the normal to the best plane
	gsl_eigen_symmv_sort(eigenValOfA, eigenVectsOfA, GSL_EIGEN_SORT_VAL_ASC);
	gsl_matrix_get_col(bestPlane, eigenVectsOfA, 0);
	//Distance from center atom to the plane is the displacement
	gsl_blas_ddot(bestPlane, centOfMass, distToPlane);

	//Cleanup
	gsl_matrix_free(eigenVectsOfA);
	gsl_vector_free(eigenValOfA);
	gsl_eigen_symmv_free(eigenA);
	gsl_matrix_free(A);
	gsl_vector_free(centOfMass);
}

void OutputStrainTensor(gsl_matrix **F, int n, string path, string fileName, string comment)
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
			ofile<<i<<" ";
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
	cout<<"Output tensor file: "<<path<<fileName<<endl;
}


