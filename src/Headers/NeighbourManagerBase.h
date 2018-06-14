/*
 * NeighbourManagerBase.h
 *
 *  Created on: 8 Dec 2016
 *      Author: rosotti
 */

#ifndef SRC_HEADERS_NEIGHBOURMANAGERBASE_H_
#define SRC_HEADERS_NEIGHBOURMANAGERBASE_H_

#include <vector>
#include <deque>
#include "TreeCell.h"
#include "Multipole.h"


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
    struct range {
      int begin, end ;
    };
	std::deque<range> tempneib;
	std::deque<range> tempperneib;
	std::deque<range> tempdirectneib;

public:
    /* Add a range of neighbours */
	template<int ndim>
    void AddNeibs(const TreeCellBase<ndim>& cell) {
      range r = {cell.ifirst, cell.ilast+1} ;
      tempneib.push_back(r);
    }
    /* Add a neighbours that we need ghosts of */
    template<int ndim>
    void AddPeriodicNeibs(const TreeCellBase<ndim>& cell) {
      range r = {cell.ifirst, cell.ilast+1} ;
      tempperneib.push_back(r);
    }
    /* Add a distant particles for gravity force */
    template<int ndim>
    void AddDirectNeibs(const TreeCellBase<ndim>& cell) {
      range r = {cell.ifirst, cell.ilast+1} ;
      tempdirectneib.push_back(r);
    }


protected:
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
	FastMultipoleForces<ndim> multipole ;
	using NeighbourManagerBase::tempneib;
    using NeighbourManagerBase::tempperneib;
    using NeighbourManagerBase::tempdirectneib;
public:
    NeighbourManagerDim() : multipole_type(monopole) { } ;

    /* Add the multipole moments of a gravity cell */
	void AddGravCell(const MultipoleMoment<ndim>& moment) {

      gravcell.push_back(moment);

	  switch (multipole_type) {
	  case monopole:
	  case quadrupole:
		break;
	  case fast_monopole:
	    multipole.AddMonopoleContribution(moment);
	    break ;
	  case fast_quadrupole:
        multipole.AddQuadrupoleContribution(moment);
        break ;
	  }
	}

	//===============================================================================================
	//  GetGravCell
	/// \brief Get the number of multipole moments stored, and return a pointer to them.
	//===============================================================================================
	int GetGravCell (MultipoleMoment<ndim>** gravcell_p) {

	  int NgravCell = gravcell.size();

	  if (NgravCell > 0)
	    *gravcell_p = &gravcell[0];
	  else
	    *gravcell_p = NULL ;

	  return NgravCell ;
	}

    //===============================================================================================
    //  ComputeFastMultipoleForces
    /// \brief Compute the fast multipole forces the given particles.
    //==============================================================================================
	template<template <int> class ParticleType>
	void ComputeFastMultipoleForces(int Nactive, ParticleType<ndim>* activepart,
	                                const ParticleTypeRegister& types) {
	  for (int j=0; j<Nactive; j++)
	    if (types[activepart[j].ptype].self_gravity)
	      multipole.ApplyForcesTaylor(activepart[j].r, activepart[j].atree, activepart[j].gpot) ;
	}


    void set_target_cell(const TreeCellBase<ndim>& cell) {
      NeighbourManagerBase::clear();
      gravcell.clear();

      multipole.set_target_cell(cell.r);
    }

    void set_multipole_type(multipole_method multipole) {
      multipole_type = multipole ;
    }

private:
    multipole_method multipole_type ;
};



#endif /* SRC_HEADERS_NEIGHBOURMANAGERBASE_H_ */
