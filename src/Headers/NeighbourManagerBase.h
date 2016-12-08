/*
 * NeighbourManagerBase.h
 *
 *  Created on: 8 Dec 2016
 *      Author: rosotti
 */

#ifndef SRC_HEADERS_NEIGHBOURMANAGERBASE_H_
#define SRC_HEADERS_NEIGHBOURMANAGERBASE_H_

#include <vector>


class NeighbourManagerBase {
protected:
	vector<int> tempneib;
	vector<int> tempperneib;
public:
	void AddPeriodicNeib(const int i) {
		tempperneib.push_back(i);
	}

	void AddNeib(const int i) {
		tempneib.push_back(i);
	}
};


#endif /* SRC_HEADERS_NEIGHBOURMANAGERBASE_H_ */
