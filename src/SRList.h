/*
 * SRList.h
 *
 *  Created on: Jun 17, 2011
 *      Author: ken
 */

#ifndef SRLIST_H_
#define SRLIST_H_

#include "Includes.h"

class SRList {
public:
	static Range** sRange;
	static int nSpaceRange;
	SRList();
	virtual ~SRList();
	Range* CreateNewSpaceRange();
	Range* FindSpaceRange(string );
	void Clean();
};

#endif /* SRLIST_H_ */
