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


//=================================================================================================
//  Class NeighbourManagerBase
/// \brief   Base class for neighbour search wrapper objects
/// \details The base classes provide the interface for storing and clearing the indices of the
///          particles found by neighbour searches. The indices should be the location of the
///          particles in the main array.
/// \author  G. Rosotti, R. A. Booth
/// \date    15/12/2016
//=================================================================================================
class NeighbourManagerBase {
protected:
	vector<int> tempneib;
	vector<int> tempperneib;
	vector<int> tempdirectneib;

public:
	/* Add a neighbour that does not need automatic ghost generation */
    void AddNeib(const int i) {
        tempneib.push_back(i);
    }
    /* Add a neighbour that we need ghosts of */
	void AddPeriodicNeib(const int i) {
		tempperneib.push_back(i);
	}
	/* Add a distant particle for gravity force */
	void AddDirectNeib(const int i) {
		tempdirectneib.push_back(i);
	}

	void clear() {
	  tempneib.clear();
	  tempperneib.clear();
	  tempdirectneib.clear();
	}
};

//=================================================================================================
//  Class NeighbourManagerDim
/// \brief   Base class for neighbour search wrapper objects, which also holds multipole data.
/// \author  G. Rosotti, R. A. Booth
/// \date    15/12/2016
//=================================================================================================
template <int ndim>
class NeighbourManagerDim : public NeighbourManagerBase {
protected:
	vector<MultipoleMoment<ndim> > gravcell;
	using NeighbourManagerBase::tempneib;
    using NeighbourManagerBase::tempperneib;
    using NeighbourManagerBase::tempdirectneib;
public:
    /* Add the multipole moments of a gravity cell */
	void AddGravCell(const MultipoleMoment<ndim>& moment) {
		gravcell.push_back(moment);
	}

	//===============================================================================================
	//  GetGravCell
	/// \brief Get the number of multipole moments stored, and return a pointer to them.
	//===============================================================================================
	int GetGravCell (MultipoleMoment<ndim>** gravcell_p) {
	  *gravcell_p = &gravcell[0];
	  return gravcell.size();
	}

    void clear() {
      NeighbourManagerBase::clear();
      gravcell.clear();
    }
};



#endif /* SRC_HEADERS_NEIGHBOURMANAGERBASE_H_ */
