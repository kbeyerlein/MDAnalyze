/*
 * VelToPhonon.cpp
 *
 *  Created on: May 14, 2010
 *      Author: Ken
 */

#include "Includes.h"
#include "Global.h"
#include "TimeSeries.h"
#include "SRList.h"

Range **SRList::sRange=0;
int SRList::nSpaceRange=0;

SRList sRList;
int nComm;
string *command=0;
TimeSeries sim;
// Design Foundations:- New Data structure which is more logical and conducive to further expansion
//					  - Use the execute as instructed philosophy, so ReadInputFile becomes program execution
//					  - Want to use as little memory as possible, means reloading data files multiple times

//NOTE: Data Structure;  TimeSeries>TimeStep>Group>VelGroup-PosGroup-Atom>pos-vel-type...

void Calc(int iComm)
{
	string c;
	stringstream temp(command[iComm]);
	temp>>c;
	if (c=="atomPos"){
		//Default parameters
		int maxKeywords=1;
		char tempBound[3];
		for (int i=0;i<3;i++){
			tempBound[i]='p'; // periodic boundary conditions
		}

		// Read in the position file identifiers
		string groupName, fileType, fileName;
		temp >> groupName>>fileType>>fileName;

		//Read in any keywords and changes to default parameters
		string keyword="";
		temp>>keyword;
		for (int i=0;i<maxKeywords&&keyword!="";i++){
			if (keyword=="boundary"){
				for (int i=0;i<3;i++){
					temp>>tempBound[i];
					if (tempBound[i]!='p'&&tempBound[i]!='m'){
						cout <<"Boundary specification "<<i<<" \'"<<tempBound[i]<<"\' not supported.";
						exit(0);
					}
				}
			}
			else{
				cout<<"Keyword: "<<keyword<<" not supported.\n";
			}
			temp>>keyword;
		}

		//Parse position file
		if (fileType=="lammpstrj"){
			sim.FindPosInLammpstrjFile(groupName, fileName, &tempBound[0]);
		}
		else{
			cout<<"Position fileType not supported: "<<fileType<<endl;
			exit(0);
		}

		//create indexed array of positions (instead of just linked list)
		sim.CreateSeriesArray();
	}
	else if (c=="atomVel"){
		string groupName, fileType, fileName;
		temp >> groupName>>fileType>>fileName;
		if (fileType=="vel"){
			sim.FindVelInLammpsVelFile(groupName, fileName);
		}
		else{
			cout<<"Velocity fileType not supported: "<<fileType<<endl;
			exit(0);
		}
		sim.CreateSeriesArray();
	}
	else if (c=="msd_t"){
		string groupName;
		int refStep, t0, tf;
		temp >> groupName>>refStep>>t0>>tf;
		sim.CalcMSD_t(groupName, refStep, t0, tf);
	}
	else if (c=="velCorrFn"){
		string groupName, timeRangeName;
		temp>>groupName>>timeRangeName;
		sim.CalcVelCorrFunc(groupName, timeRangeName);
	}
	else if (c=="disp"){
		string type;
		temp>>type;
		if (type=="specific"){
			string groupName, remaining;
			int t0, tf;
			temp>>groupName>>t0>>tf;
			getline(temp, remaining);
			Timestep *t_0=sim.FindTimestep(t0), *t_f=sim.FindTimestep(tf);
			if (t_0!=0&&t_f!=0){
				t_0->InputPos("all");
				t_f->InputPos("all");
			}
			else{
				exit(0);
			}
			if (remaining.find("scale")!=remaining.npos){
				sim.CalcDisplace(&sim.userDisp, "timeStep", groupName, true, t0, tf);
			}
			else{
				sim.CalcDisplace(&sim.userDisp, "timeStep", groupName, false, t0, tf);
			}
			t_0->Clean();
			t_f->Clean();
			if (remaining.find("pMSD")!=remaining.npos){
				sim.userDisp->CalcPMSD();
			}
			if (remaining.find("radMSD")!=remaining.npos){
				Timestep *t=sim.FindTimestep(sim.userDisp->initial);
				if (t!=0){
					t->InputPos("all");
				}
				else{
					exit(0);
				}
				t->all.CalcCenter();
				sim.userDisp->CalcRadMSD(t->all.p->center);
				t->Clean();
			}
			if (remaining.find("radDisp")!=remaining.npos){
				Timestep *t=sim.FindTimestep(sim.userDisp->initial);
				if (t!=0){
					t->InputPos("all");
				}
				else{
					exit(0);
				}
				t->all.CalcCenter();
				sim.userDisp->CalcRadDisp(t->all.p->center);
				t->Clean();
			}
		}
		else if (type=="static"){
			string groupName, remaining;
			int ref, t0, tf;
			temp>>groupName>>ref>>t0>>tf;
			getline(temp, remaining);
			if (remaining.find("scale")!=remaining.npos){
				cout<<"Scaled displacement currently not supported.\n";
				exit(0);
				//sim.CalcStaticDisp(groupName, true, ref, t0, tf);
			}
			else{
				sim.CalcStaticDisp(groupName, false, ref, t0, tf);
			}
			if (remaining.find("pMSD")!=remaining.npos){
				sim.staticDisp->CalcPMSD();
			}
			if (remaining.find("radMSD")!=remaining.npos){
				Timestep *t=sim.FindTimestep(sim.staticDisp->initial);
				if(t!=0){
					t->InputPos("all");
				}
				else{
					exit(0);
				}
				t->all.CalcCenter();
				sim.staticDisp->CalcRadMSD(t->all.p->center);
				t->Clean();
			}
			if (remaining.find("radDisp")!=remaining.npos){
				Timestep *t=sim.FindTimestep(sim.staticDisp->initial);
				if (t==0){
					exit(0);
				}
				t->InputPos("all");
				t->all.CalcCenter();
				sim.staticDisp->CalcRadDisp(t->all.p->center);
				t->Clean();
			}
		}
		else if (type=="dynamic"){
			bool pMSD=false, radMSD=false, radDisp=false;
			string groupName, remaining;
			int t0, tf, dt;
			temp>>groupName>>t0>>tf>>dt;
			getline(temp, remaining);
			if (remaining.find("scale")!=remaining.npos){
				sim.CalcDynamicDisp(groupName, true, t0, tf, dt);
			}
			else{
				sim.CalcDynamicDisp(groupName, false, t0, tf, dt);
			}
			if (remaining.find("pMSD")!=remaining.npos){
				pMSD=true;
			}
			if (remaining.find("radMSD")!=remaining.npos){
				radMSD=true;
			}
			if (remaining.find("radDisp")!=remaining.npos){
				radDisp=true;
			}
			for (int i=0;i<sim.nDynDisp;i++){
				if (pMSD){
					sim.dynDisp[i]->CalcPMSD();
				}
				if (radMSD){
					sim.avg->CalcCenter();
					sim.dynDisp[i]->CalcRadMSD(sim.avg->p->center);
				}
				if (radDisp){
					sim.avg->CalcCenter();
					sim.dynDisp[i]->CalcRadDisp(sim.avg->p->center);
				}
			}
		}
		else{
			cout<<"Displacement type is not supported: "<<type<<endl;
		}
	}
	else if (c=="atomStrain"){
		string groupName, weightType;
		int t0, tf, dt;
		double Rc;
		temp>>groupName>>t0>>tf>>dt>>Rc>>weightType;
		sim.CalcAtomicStrainTensors(groupName, t0, tf, dt, Rc, weightType);
	}
	else if (c=="avgPos"){
		string groupName, tRangeName;
		temp>>tRangeName>>groupName;
		Range *tR=sim.FindTimeRange(tRangeName);
		sim.CalcAvgPos(groupName, tR->x0, tR->xf);
	}
	else{
		cout<<"Command: "<<c<<" is not supported.\n";
	}

}
void Define(int iComm)
{
	string d;
	stringstream temp(command[iComm]);
	temp>>d>>d;
	if (d=="neighRange"){
		string groupName1, groupName2;
		int ref;
		double r0,rc;
		temp>>ref>>groupName1>>groupName2>>r0>>rc;

		Group *g1, *g2;
		Timestep *t=0;
		if (ref==-1){
			if (!sim.avg){
				cout<<"Cannot calculate Avg Atom Distance in Range for Avg because it has not been calculated.\n";
				exit(0);
			}
			if(groupName1=="all"&&groupName2=="all"){
				g1=sim.avg;
				g2=sim.avg;
			}
			else {
				cout<<"Groups of average besides 'all' are no currently supported.\n";
				exit(0);
			}
		}
		else{
			t=sim.FindTimestep(ref);
			if (t==0){
				exit(0);
			}
			g1=t->FindGroup(groupName1);
			g2=t->FindGroup(groupName2);
		}
		g1->neighR0=r0;
		g1->neighRc=rc;
		g1->neighG=g2;
	}
	else if(d=="timeRange"){
		Range *tR=sim.CreateNewTimeRange();
		temp>>tR->name>>tR->x0>>tR->xf>>tR->dx;
	}
	else if (d=="spaceRange"){
		SRList list;
		Range *sR=list.CreateNewSpaceRange();
		temp>>sR->name>>sR->x0>>sR->xf>>sR->dx;
	}
	else if (d=="newTimestep"){
		string type;
		temp>>type;
		if (type=="move"){
			string dataType;
			int tOld, tNew;
			string nameOld, nameNew;
			temp>>dataType>>tOld>>nameOld>>tNew>>nameNew;
			if (dataType=="pos"){
				sim.MovePosToNewTime(tOld, nameOld, tNew, nameNew);
			}
			else{
				cout<<"newTimestep move data type: "<<dataType<< " is not supported.\n";
			}
		}
		else{
			cout<<"newTimestep creation type: "<<type<< " is not supported.\n";
		}
	}
	else if (d=="posInputAction"){
		string tRName, action;
		temp>>tRName;
		getline(temp,action);
		Range *timeRange=sim.FindTimeRange(tRName);
		int tMin=timeRange->x0;
		for (int i=0;i<sim.nTimeSteps;i++){
			if(sim.step[i]->step>=tMin && sim.step[i]->step<=timeRange->xf){
				sim.step[i]->CreatePosInputAction(action);
				tMin+=timeRange->dx;
			}
		}
	}
	else if (d=="group"){
		string tRName, groupName, groupCommand;
		temp>>groupName>>tRName;
		getline(temp,groupCommand);
		sim.DefineUserGroup(tRName, groupName, groupCommand);
	}
	else{
		cout<<"Command: "<<d<<" is not supported.\n";
	}
}
void Output(int iComm)
{
	string o;
	stringstream temp(command[iComm]);
	temp>>o>>o;
	if(o=="msd_t"){
		string path, name;
		temp>>path>> name;
		if (sim.msd_t){
			OutputDistrib(sim.msd_t, path, name, false, false, "");
		}
	}
	else if (o=="velCorrFn"){
		string path, name;
		temp >> path>>name;
		if (sim.velCorrFn){
			OutputDistrib(sim.velCorrFn, path, name, false, true, "");
		}
	}
	else if (o=="disp"){
		string path, dist, name, type;
		temp>>type>>dist>>path>>name;
		if (type=="specific"){
			if (sim.userDisp){
				if (dist=="pMSD"){
					if (sim.userDisp->pMSD){
						OutputDistrib(sim.userDisp->pMSD, path, name, false, true, "");
					}
				}
				else if (dist=="radMSD"){
					if (sim.userDisp->radMSD){
						OutputDistrib(sim.userDisp->radMSD, path, name, true, false, "");
					}
				}
				else if (dist=="radDisp"){
					if (sim.userDisp->radDisp){
						OutputDistrib(sim.userDisp->radDisp, path, name, true, false, "");
					}
				}
			}
		}
		else if (type=="static"){
			if (sim.staticDisp){
				if (dist=="pMSD"){
					if (sim.staticDisp->pMSD){
						OutputDistrib(sim.staticDisp->pMSD, path, name, false, true, "");
					}
				}
				else if (dist=="radMSD"){
					if (sim.staticDisp->radMSD){
						OutputDistrib(sim.staticDisp->radMSD, path, name, true, false, "");
					}
				}
				else if (dist=="radDisp"){
					if (sim.staticDisp->radDisp){
						OutputDistrib(sim.staticDisp->radDisp, path, name, true, false, "");
					}
				}
			}
		}
		else if (type=="dynamic"){
			for (int i=0;i<sim.nDynDisp;i++){
				if (sim.dynDisp[i]){
					stringstream comment;
					comment<<"Distribution from displacement calculated between ";
					if (sim.dynDisp[i]->initial==-1){
						comment<<"Avg. Pos";
					}
					else{
						comment<<"Step "<<sim.dynDisp[i]->initial;
					}
					comment <<" and ";
					if (sim.dynDisp[i]->final==-1){
						comment<<"Avg. Pos";
					}
					else{
						comment<<"Step "<<sim.dynDisp[i]->final;
					}
					if (dist=="pMSD"){
						if (sim.dynDisp[i]->pMSD){
							OutputDistrib(sim.dynDisp[i]->pMSD, path, name, false, true, comment.str());
						}
					}
					else if (dist=="radMSD"){
						if (sim.dynDisp[i]->radMSD){
							OutputDistrib(sim.dynDisp[i]->radMSD, path, name, true, false, comment.str());
						}
					}
					else if (dist=="radDisp"){
						if (sim.dynDisp[i]->radDisp){
							OutputDistrib(sim.dynDisp[i]->radDisp, path, name, true, false, comment.str());
						}
					}
				}
			}
		}
	}
	else if(o=="atomStrain"){
		string type, path, name;
		temp>>type>>path>>name;
		if (type=="deformGrad"){
			for (int i=0;i<sim.nStrTen;i++){
				stringstream t;
				t<<"Atomic Deformation Gradient between timesteps: "<<sim.F[i].t0<<" and "<<sim.F[i].tf;
				sim.OutputStrainTensor(&sim.F[i], path, name, t.str());
			}
		}
		else if(type=="Lagrange-Green"){
			for (int i=0;i<sim.nStrTen;i++){
				AtomStrain E;
				sim.CalcLGStrainTensor(i,&E);
				stringstream t;
				t<<"Atomic Lagrange-Green strain between timesteps: "<<sim.F[i].t0<<" and "<<sim.F[i].tf;
				sim.OutputStrainTensor(&E, path, name, t.str());
				sim.CleanAtomTensor(&E);
			}
		}
		else if(type=="Euler-Almansi"){
			for (int i=0;i<sim.nStrTen;i++){
				AtomStrain e;
				sim.CalcEAStrainTensor(i,&e);
				stringstream t;
				t<<"Atomic Euler-Almansi strain between timesteps: "<<sim.F[i].t0<<" and "<<sim.F[i].tf;
				sim.OutputStrainTensor(&e, path, name, t.str());
				sim.CleanAtomTensor(&e);
			}
		}
	}
	else if (o=="avgStrain"){
		string type, path, name;
		temp>>type>>path>>name;
		gsl_matrix *AvgF=gsl_matrix_calloc(3,3);
		int t0=0;
		for (int i=0;i<sim.nStrTen;i++){
			stringstream t;
			if (i==0){
				t0=sim.F[0].t0;
			}
			gsl_matrix *tempAvg=gsl_matrix_calloc(3,3);
			sim.CalcAvgAtomStrain(i,tempAvg);
			if (i==0){
				gsl_matrix_add(AvgF,tempAvg);
			}
			else{
				//compute F_tot=F_n*F_n-1*...F_1
				gsl_matrix *tempTot=gsl_matrix_alloc(3,3);
				gsl_matrix_memcpy(tempTot, AvgF);
				gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, tempAvg, tempTot, 0.0, AvgF);
				gsl_matrix_free(tempTot);
			}
			if (type=="defGrad"){
				t<<"Avg. Atomic Deformation Gradient between timesteps: "<<sim.F[i].t0<<" and "<<sim.F[i].tf;
				OutputStrainTensor(&AvgF, 1, path, name, t.str());
			}
			else if(type=="rot"||type=="U"||type=="V"){
				gsl_matrix *R=gsl_matrix_calloc(3,3);
				gsl_matrix *U=gsl_matrix_calloc(3,3);
				gsl_matrix *V=gsl_matrix_calloc(3,3);
				sim.DecompDefGradTensor(AvgF, R, U, V);
				if (type=="rot"){
					t<<"Rigid body rotation between timesteps: "<<sim.F[i].t0<<" and "<<sim.F[i].tf;
					OutputStrainTensor(&R, 1, path, name, t.str());
				}
				if (type=="U"){
					t<<"Avg. Right Stretch Tensor between timesteps: "<<sim.F[i].t0<<" and "<<sim.F[i].tf;
					OutputStrainTensor(&U, 1, path, name, t.str());
				}
				if (type=="V"){
					t<<"Avg. Left Stretch Tensor between timesteps: "<<sim.F[i].t0<<" and "<<sim.F[i].tf;
					OutputStrainTensor(&V, 1, path, name, t.str());
				}
				gsl_matrix_free(R);
				gsl_matrix_free(U);
				gsl_matrix_free(V);

			}
			gsl_matrix_free(tempAvg);
		}
		gsl_matrix_free(AvgF);
	}
	else if (o=="atomVolChange"){
		string path, name;
		temp>>path>>name;
		for (int i=0;i<sim.nStrTen;i++){
			double *dV=0;
			try{
				dV = new double [sim.F[i].nA];
			}
			catch(exception &e){

			}
			for (int j=0;j<sim.F[i].nA;j++){
				dV[j]=CalcJacobian_gsl(sim.F[i].s[j]);
			}
			Timestep *curTime=sim.FindTimestep(sim.F[i].t0);
			if (curTime==0){
				exit(0);
			}
			curTime->InputPos("all");
			sim.F[i].g_i->OutputVolChangeMap(dV, path, name, sim.F[i].t0, i+1);
			curTime->Clean();
			delete [] dV;
		}
	}
	else if(o=="avgDistRDF"){
		int time;
		string g1,g2, sRangeName, path, name;
		Timestep *curTime=0;
		Group *curGroup, *neighGroup;
		Range *sR;
		Distrib RDF;
		temp>>time>>sRangeName>>g1>>g2>>path>>name;
		if (time!=-1){
			curTime=sim.FindTimestep(time);
			if (curTime==0){
				exit(0);
			}
			curTime->InputPos("all");
			curGroup=curTime->FindGroup(g1);
			neighGroup=curTime->FindGroup(g2);
		}
		else{
			if (g1=="all"&&g2=="all"){
				curGroup=sim.avg;
				neighGroup=sim.avg;
			}
			else{
				cout<<"Groups of average besides 'all' are no currently supported.\n";
				exit(0);
			}
		}
		sR=sRList.FindSpaceRange(sRangeName);
		AllocDistrib(&RDF, "RDF", "Dist.", "Num", sR->x0, sR->xf, sR->dx);
		curGroup->BuildNeighsInRange_fast(neighGroup, sR->x0, sR->xf);
		curGroup->CalcRDF(&RDF);
		OutputDistrib(&RDF, path, name, false, true, "");
		CleanDistrib(&RDF);
		if (curTime){
			curTime->Clean();
		}
	}
	else if(o=="avgNeighDist_R"){
		int time;
		double dr;
		string g1, g2, sRangeName, path, name;
		Timestep *curTime=0;
		Group *curGroup, *neighGroup;
		Range *sR;
		Distrib NN_R;
		temp>>time>>sRangeName>>g1>>g2>>dr>>path>>name;

		if (time!=-1){
			curTime=sim.FindTimestep(time);
			if (curTime==0){
				exit(0);
			}
			curTime->InputPos("all");
			curGroup=curTime->FindGroup(g1);
			neighGroup=curTime->FindGroup(g2);
		}
		else{
			if (g1=="all"&&g2=="all"){
				curGroup=sim.avg;
				neighGroup=sim.avg;
			}
			else{
				cout<<"Groups of average besides 'all' are no currently supported.\n";
				exit(0);
			}
		}
		sR=sRList.FindSpaceRange(sRangeName);
		curGroup->BuildNeighsInRange_fast(neighGroup, sR->x0, sR->xf);
		curGroup->CalcCenter();
		curGroup->CalcPosMags(curGroup->p->center);
		AllocDistrib(&NN_R, "NeighDist(R)", "R", "Dist", 0, curGroup->GetMaxR(), dr);
		curGroup->CalcAvgNeighDistVsR(&NN_R);
		OutputDistrib(&NN_R, path, name, true, true, "");
		CleanDistrib(&NN_R);
		if (curTime){
			curTime->Clean();
		}
	}
	else if (o=="avgNeighDist_t"){
		string tR,sR,g1,g2,path, name;
		Distrib NN_t;
		temp>>tR>>sR>>g1>>g2>>path>>name;
		sim.CalcNN_t(&NN_t, tR, sR, g1,g2);
		OutputDistrib(&NN_t, path, name, false, true, "");
		CleanDistrib(&NN_t);
	}
	else if (o=="numNN"){
		int time;
		double maxNN;
		string g1, g2, sRangeName, path, name;
		Timestep *curTime=0;
		Group *curGroup, *neighGroup;
		Range *sR;
		Distrib nNN;
		temp>>time>>sRangeName>>g1>>g2>>maxNN>>path>>name;

		if (time!=-1){
			curTime=sim.FindTimestep(time);
			if (curTime==0){
				exit(0);
			}
			curTime->InputPos("all");
			curGroup=curTime->FindGroup(g1);
			neighGroup=curTime->FindGroup(g2);
		}
		else{
			if (g1=="all"&&g2=="all"){
				curGroup=sim.avg;
				neighGroup=sim.avg;
			}
			else{
				cout<<"Groups of average besides 'all' are not currently supported.\n";
				exit(0);
			}
		}
		sR=sRList.FindSpaceRange(sRangeName);
		//curGroup->BuildNeighsInRange(neighGroup, sR->x0, sR->xf);
		curGroup->BuildNeighsInRange_fast(neighGroup, sR->x0, sR->xf);
		AllocDistrib(&nNN, "numNeigh", "nNN", "nAtoms", 0, maxNN, 1.0);
		curGroup->CalcNumNNDist(&nNN, maxNN);
		OutputDistrib(&nNN, path, name, false, true, "");
		CleanDistrib(&nNN);
		if (curTime){
			curTime->Clean();
		}
	}
	else if (o=="avgPos"){
		//the average position must already be calculated previously.
		string format, path, name;
		temp>>format>>path>>name;
		sim.avg->OutputPosition(format,path,name);

	}
	else if (o=="avgRadius"){
		string timeRangeName, groupName, path, fileName;
		temp>>timeRangeName>>groupName>>path>> fileName;
		sim.CalcRadius(timeRangeName, groupName, path, fileName);
	}
	else if (o=="densityVsR"){
		int time;
		double dr;
		string g1, path, name;
		Timestep *curTime=0;
		Group *curGroup;
		Distrib density;
		gsl_vector *center=0;
		temp>>time>>g1>>dr>>path>>name;
		if (time!=-1){
			curTime=sim.FindTimestep(time);
			if (curTime==0){
				exit(0);
			}
			curTime->InputPos("all");
			curGroup=curTime->FindGroup(g1);
			curTime->all.CalcCenter();
			center=curTime->all.p->center;
		}
		else{
			if (g1=="all"){
				curGroup=sim.avg;
				sim.avg->CalcCenter();
				center=sim.avg->p->center;

			}
			else{
				cout<<"Groups of average besides 'all' are not currently supported.\n";
				exit(0);
			}
		}
		curGroup->CalcPosMags(center);
		double maxR=curGroup->GetMaxR()+dr;
		AllocDistrib(&density, "DensityVsR", "R", "density", 0, maxR, dr);
		curGroup->CalcDensityVsR(&density);
		OutputDistrib(&density, path, name, true, true, "");
		CleanDistrib(&density);
		if (curTime){
			curTime->Clean();
		}

	}
	else{
		cout<<"Output of calc not supported: "<<o<<endl;
		exit(0);
	}
}
void ReadInputFile(string dir, string fileName)
{
	string temp;
	string file=dir+"/"+fileName;
	ifstream input;
	if (nComm==0){
		input.open(file.c_str(), ios::in);
		if(input.is_open()){
			while(!input.eof()){
				getline(input, temp);
				nComm++;
			}
			input.close();
			ReadInputFile(dir, fileName);
		}
		else{
			cout<<"Cannot open Input file: "<<file<<endl;
			exit(0);
		}
	}
	else{
		try{
			if (command){
				delete [] command;
			}
			command=new string [nComm];
		}
		catch(exception &e){
			cout<<e.what()<<endl;
		}
		input.open(file.c_str(), ios::in);
		if(input.is_open()){
			int count=0;
			cout<<"Reading from input file: "<<file<<endl;
			while(!input.eof()){
				getline(input, command[count]);
				count++;
			}
			input.close();
		}
	}
}


int main (int argc, char **argv)
{
	string infile="", path="";
	for (int i = 0; i < argc; i++){
		if (string(argv[i])=="-f"){
			infile = argv[i+1];
		}
		if (string(argv[i])=="-p"){
			path=argv[i+1];
		}
	}
	cout<<"MDAnalyze v1.0  -Nov 2010-  Created by: Kenneth R. Beyerlein\n";
	// If not input file is specified
	if (infile == "")
	{
		cout << "No input file specified...Exiting\n";
		return 0;
	}
	// If no directory specified
	if (path==""){
		//path="C:/Users/Ken/Documents/Input/VelToPhonon";
		path=".";
		cout<<"Default input file path assumed: " <<path<<endl;
	}
	// Reading from the input file and Initializing parameters
	ReadInputFile(path, infile);
	for(int i=0;i<nComm;i++){
		string tempComm;
		stringstream temp (command[i]);
		temp>>tempComm;
		if (tempComm=="name"){
			temp>>sim.name;
		}
		else if(tempComm=="output"){
			Output(i);
		}
		else if (tempComm=="define"){
			Define(i);
		}
		else if (tempComm==""||tempComm=="#"){
		}
		else{
			Calc(i);
		}
	}
	if (command){
		delete [] command;
	}
	return 0;
}
