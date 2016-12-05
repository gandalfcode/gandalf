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

class NeighbourManagerBase {
protected:
	vector<int> tempneib;
public:
	void add_tempneib(const int i) {
		tempneib.push_back(i);
	}
};


template <int ndim, template<int> class ParticleType>
class NeighbourManager : public NeighbourManagerBase {
private:
	vector<int> neiblist;
	vector<ParticleType<ndim> > neibpart;

	vector<int> sphlist;

	vector<FLOAT> dr;
	vector<FLOAT> drmag;

	DomainBox<ndim> _domain;
public:
	using NeighbourManagerBase::tempneib;

	void clear() {
		neiblist.clear();
		tempneib.clear();
		neibpart.clear();
	}

	void EndSearch(const TreeCellBase<ndim> &cell, const ParticleType<ndim>* partdata, const DomainBox<ndim>& domain,
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
		  for (int ii=0; ii < tempneib.size(); ii++) {
		    const int i = tempneib[ii] ;


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


	}

	int GetNeibList(const ParticleType<ndim>& p,const Typemask& hydromask, int** neiblist_p, ParticleType<ndim>** neibpart_p,
			FLOAT** drmag_p, FLOAT** dr_p, int* levelneib) {


	    FLOAT rp[ndim];                                // Local copy of particle position
	    FLOAT draux[ndim];                             // Aux. relative position vector
	    FLOAT drsqd;                                   // Distance squared


        for (int k=0; k<ndim; k++) rp[k] = p.r[k];

        const FLOAT hrangesqdi = p.hrangesqd;

        const GhostNeighbourFinder<ndim> GhostFinder(_domain);


        // Validate that gather neighbour list is correct
#if defined(VERIFY_ALL)
        if (neibcheck) this->CheckValidNeighbourList(i, sph->Nhydro + sph->NPeriodicGhost,
                                                     Nneib, neiblist, sphdata, "all");
#endif

        dr.clear();
        drmag.clear();
        sphlist.clear();

        // Compute distances and the inverse between the current particle and all neighbours here,
        // for both gather and inactive scatter neibs.  Only consider particles with j > i to
        // compute pair forces once unless particle j is inactive.
        //-----------------------------------------------------------------------------------------
        for (int jj=0; jj<neibpart.size(); jj++) {

          // Skip non-hydro particles and the current active particle.
          if (hydromask[neibpart[jj].ptype] == false) continue;

          for (int k=0; k<ndim; k++) draux[k] = neibpart[jj].r[k] - rp[k];
          GhostFinder.NearestPeriodicVector(draux);
          drsqd = DotProduct(draux, draux, ndim) + small_number;

          if (drsqd <= small_number) continue ;

          // Compute relative position and distance quantities for pair
          if (drsqd <= hrangesqdi || drsqd <= neibpart[jj].hrangesqd) {
        	const FLOAT drmagi = sqrt(drsqd);
            drmag.push_back(drmagi);
            const FLOAT invdrmag = (FLOAT) 1.0/drmagi;
            for (int k=0; k<ndim; k++) dr.push_back(draux[k]*invdrmag);
            levelneib[neiblist[jj]] = max(levelneib[neiblist[jj]],p.level);
            sphlist.push_back(jj);
          }
        }

        *neiblist_p=&sphlist[0];
        *neibpart_p=&neibpart[0];
        *drmag_p=&drmag[0];
        *dr_p=&dr[0];

        assert(drmag.size()==sphlist.size());

        return sphlist.size();

	}

};


#endif /* SRC_HEADERS_NEIGHBOURMANAGER_H_ */
