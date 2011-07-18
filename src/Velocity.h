/*
 * Velocity.h
 *
 *  Created on: May 14, 2010
 *      Author: Ken
 */

#ifndef VELOCITY_H_
#define VELOCITY_H_
#include "Includes.h"
#include "Global.h"
class Velocity {
public:
	Velocity();
	virtual ~Velocity();

	int dataPos;
	string file, fileType;
	Distrib *pVel;

};

#endif /* VELOCITY_H_ */
