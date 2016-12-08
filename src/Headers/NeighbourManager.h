/*
 * NeighbourManager.h
 *
 *  Created on: 5 Dec 2016
 *      Author: rosotti
 */

#ifndef SRC_HEADERS_NEIGHBOURMANAGER_H_
#define SRC_HEADERS_NEIGHBOURMANAGER_H_


#include <vector>
using namespace std;
#include "TreeCell.h"
#include "GhostNeighbours.hpp"
#include "NeighbourManagerBase.h"



template <int ndim, class ParticleType>
class NeighbourManager : public NeighbourManagerBase {
private:
	int _NPeriodicGhosts;
	vector<int> neiblist;
	vector<ParticleType > neibpart;

	vector<int> culled_neiblist;

	vector<FLOAT> dr;
	vector<FLOAT> drmag;

	DomainBox<ndim> _domain;
public:
	using NeighbourManagerBase::tempneib;

	void clear() {
		neiblist.clear();
		tempneib.clear();
		tempperneib.clear();
		neibpart.clear();
	}

	template<class InParticleType>
	void EndSearch(const TreeCellBase<ndim> &cell, const InParticleType* partdata, const DomainBox<ndim>& domain,
			const FLOAT kernrange) {

		  assert(partdata != NULL);

		  _domain = domain;

		  FLOAT dr[ndim];                      // Relative position vector
		  FLOAT drsqd;                         // Distance squared
		  FLOAT rc[ndim];                      // Position of cell
		  const FLOAT hrangemaxsqd = pow(cell.rmax + kernrange*cell.hmax,2);
		  const FLOAT rmax = cell.rmax;
		  for (int k=0; k<ndim; k++) rc[k] = cell.rcell[k];


		  const GhostNeighbourFinder<ndim> GhostFinder(domain, cell) ;
		  int MaxGhosts = GhostFinder.MaxNumGhosts() ;

		  int Nneib = 0;

		  // Now load the particles
		  // First the ones that need ghosts to be created on the fly
		  for (int ii=0; ii < tempperneib.size(); ii++) {
		    const int i = tempperneib[ii] ;
			if (partdata[i].flags.is_dead()) continue;


		    const int NumGhosts = GhostFinder.ConstructGhostsScatterGather(partdata[i], neibpart);

		    int Nmax = NumGhosts + Nneib;
		    while (Nneib < Nmax) {
		      //neibpart[Nneib].iorig = i;
		      for (int k=0; k<ndim; k++) dr[k] = neibpart[Nneib].r[k] - rc[k];
		      drsqd = DotProduct(dr, dr, ndim);
		      FLOAT h2 = rmax + kernrange*neibpart[Nneib].h;
		      if (drsqd < hrangemaxsqd || drsqd < h2*h2) {
		        neiblist.push_back(i);
		        Nneib++ ;
		      }
		      else if (Nmax > Nneib) {
		        Nmax-- ;
		        if (Nmax > Nneib)
		          neibpart[Nneib] = neibpart[Nmax];
		          neibpart.resize(neibpart.size()-1);
		      }
		    }// Loop over Ghosts
		  }

		  assert(Nneib==neiblist.size());
		  assert (Nneib==neibpart.size());

		  _NPeriodicGhosts = Nneib;


		  for (int ii=0; ii<tempneib.size(); ii++) {
			  const int i = tempneib[ii];
			  if (partdata[i].flags.is_dead()) continue;

			  for (int k=0; k<ndim; k++) dr[k] = partdata[i].r[k] - rc[k];
			    drsqd = DotProduct(dr,dr,ndim);
			    if (drsqd < hrangemaxsqd || drsqd <
			        (rmax + kernrange*partdata[i].h)*(rmax + kernrange*partdata[i].h)) {
			      neibpart.push_back(partdata[i]);
			      neiblist.push_back(i);
			      Nneib++;
			    }
		  }

		  assert(Nneib==neiblist.size());
		  assert(neiblist.size()==neibpart.size());


	}

	int GetNumAllNeib() {
		return neiblist.size();
	}

	std::pair<int,ParticleType*> GetNeibI(int i) {
		return make_pair(neiblist[i],&neibpart[i]);
	}

	template<class InParticleType>
	int GetParticleNeib(const InParticleType& p,const Typemask& hydromask, int** neiblist_p,
			ParticleType** neibpart_p, const bool do_pair_once) {


	    FLOAT rp[ndim];                                // Local copy of particle position
	    FLOAT draux[ndim];                             // Aux. relative position vector
	    FLOAT drsqd;                                   // Distance squared


        for (int k=0; k<ndim; k++) rp[k] = p.r[k];

        const FLOAT hrangesqdi = p.hrangesqd;

        const GhostNeighbourFinder<ndim> GhostFinder(_domain);

        culled_neiblist.clear();

        // Compute distances and the inverse between the current particle and all neighbours here,
        // for both gather and inactive scatter neibs.  Only consider particles with j > i to
        // compute pair forces once unless particle j is inactive.
        //-----------------------------------------------------------------------------------------
        for (int jj=0; jj<neibpart.size(); jj++) {


			// Skip if (i) neighbour particle type does not interact hydrodynamically with particle,
			// (ii) neighbour is a dead (e.g. accreted) particle. Additionally, if do_pair_once is true,
        	/// (iii) same i.d. as current active particle, (iv) neighbour is on lower timestep level (i.e.
        	// timestep is shorter), or (v) neighbour is on same level as current particle but has larger id. value
			// (to only calculate each pair once).
        	if (hydromask[neibpart[jj].ptype] == false) continue;

        	if (p.iorig==neibpart[jj].iorig) continue;

        	if (do_pair_once) {

			  const bool need_interaction =
								  (neibpart[jj].flags.is_mirror()) ||
								  (p.level > neibpart[jj].level) ||
								  (p.level == neibpart[jj].level &&
								   p.iorig < neibpart[jj].iorig) ;

			  if (not need_interaction) continue;
        	}

          for (int k=0; k<ndim; k++) draux[k] = neibpart[jj].r[k] - rp[k];
          if (jj<_NPeriodicGhosts) {
			  GhostFinder.NearestPeriodicVector(draux);
          }
		  drsqd = DotProduct(draux,draux,ndim)+small_number;

//          if (drsqd <= small_number) continue ;

          // Compute relative position and distance quantities for pair
          if (drsqd <= hrangesqdi || drsqd <= neibpart[jj].hrangesqd) {
            culled_neiblist.push_back(jj);
            if (jj<_NPeriodicGhosts) for (int k=0; k<ndim; k++) neibpart[jj].r[k] = draux[k]+rp[k];
          }
        }

        *neiblist_p=&culled_neiblist[0];
        *neibpart_p=&neibpart[0];



        return culled_neiblist.size();

	}

};


#endif /* SRC_HEADERS_NEIGHBOURMANAGER_H_ */
