/*
 * NeighbourManagerBase.h
 *
 *  Created on: 8 Dec 2016
 *      Author: rosotti
 */

#ifndef SRC_HEADERS_NEIGHBOURMANAGERBASE_H_
#define SRC_HEADERS_NEIGHBOURMANAGERBASE_H_

#include <vector>
#include "TreeCell.h"

class NeighbourManagerBase {
protected:
	vector<int> tempneib;
	vector<int> tempperneib;
	vector<int> tempdirectneib;

public:
	void AddPeriodicNeib(const int i) {
		tempperneib.push_back(i);
	}

	void AddNeib(const int i) {
		tempneib.push_back(i);
	}

	void AddDirectNeib(const int i) {
		tempdirectneib.push_back(i);
	}
};

template <int ndim>
class NeighbourManagerDim : public NeighbourManagerBase {
protected:
	vector<MultipoleMoment<ndim> > gravcell;
	using NeighbourManagerBase::tempneib;
    using NeighbourManagerBase::tempperneib;
    using NeighbourManagerBase::tempdirectneib;
public:
	void AddGravCell(const MultipoleMoment<ndim>& moment) {
		gravcell.push_back(moment);
	}
};



#endif /* SRC_HEADERS_NEIGHBOURMANAGERBASE_H_ */
