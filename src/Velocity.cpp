/*
 * Velocity.cpp
 *
 *  Created on: May 14, 2010
 *      Author: Ken
 */

#include "Velocity.h"

Velocity::Velocity() {
	pVel=0;
	file="";
	dataPos=0;

}

Velocity::~Velocity() {
	if (pVel){
		CleanDistrib(pVel);
		delete pVel;
	}
}




