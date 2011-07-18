/*
 * SRList.cpp
 *
 *  Created on: Jun 17, 2011
 *      Author: ken
 */

#include "SRList.h"

SRList::SRList() {
	// No member class variables

}

SRList::~SRList() {
	// No member class variables, they are all static.
	// Cleaning of SRList is taken care of in TimeSeries destructor.
}
Range* SRList::CreateNewSpaceRange()
{
	Range **temp=0;
	if (sRange){
		if (nSpaceRange>1){
			temp=new Range* [nSpaceRange];
			for(int i=0;i<nSpaceRange;i++){
				temp[i]=sRange[i];
			}
			delete [] sRange;
		}
		else{
			temp=new Range*;
			*temp=*sRange;
			delete sRange;
		}
	}
	nSpaceRange++;
	if (nSpaceRange>1){
		sRange=new Range* [nSpaceRange];
		for (int i=0;i<nSpaceRange-1;i++){
			sRange[i]=temp[i];
		}
		if(nSpaceRange==2){
			delete temp;
		}
		else{
			delete [] temp;
		}
		sRange[nSpaceRange-1]=new Range;
		return sRange[nSpaceRange-1];
	}
	else{
		sRange=new Range *;
		*sRange=new Range;
		return *sRange;
	}
}
Range* SRList::FindSpaceRange(string id)
{
	for (int i=0;i<nSpaceRange;i++){
		if (sRange[i]->name==id){
			return sRange[i];
		}
	}
	cout<<"\nSpaceRange '"<<id<<"' has not been defined.\n";
	return 0;
}
void SRList::Clean()
{
	if (sRange){
		for (int i=0;i<nSpaceRange;i++){
			delete sRange[i];
		}
		if (nSpaceRange>1){
			delete [] sRange;
		}
		else{
			delete sRange;
		}
	}
}
