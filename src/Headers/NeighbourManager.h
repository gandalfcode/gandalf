/*
 * NeighbourManager.h
 *
 *  Created on: 5 Dec 2016
 *      Author: rosotti
 */

#ifndef NEIGHBOURMANAGER_H_
#define NEIGHBOURMANAGER_H_


#include <vector>
using namespace std;
#include "TreeCell.h"
#include "GhostNeighbours.hpp"
#include "NeighbourManagerBase.h"
#include "Particle.h"

struct ListLength {
  int Nhydro;
  int Ndirect;
  int Ngrav;
};

template <int ndim, class ParticleType>
class NeighbourManager : public NeighbourManagerDim<ndim> {
private:
  int _NPeriodicGhosts;
  int _NCellDirectNeib;
  vector<int> neiblist;
  vector<int> directlist;
  vector<int> neib_idx ;
  vector<ParticleType > neibdata;

  vector<int> culled_neiblist;
  vector<int> gravlist;

  vector<FLOAT> dr;
  vector<FLOAT> drmag;

  const DomainBox<ndim>* _domain;
  ParticleTypeRegister _types ;
  double _kernrange ;

public:
  using NeighbourManagerDim<ndim>::tempneib;
  using NeighbourManagerDim<ndim>::tempperneib;
  using NeighbourManagerDim<ndim>::tempdirectneib;
  using NeighbourManagerDim<ndim>::gravcell;

  NeighbourManager(const Hydrodynamics<ndim>* hydro, const DomainBox<ndim>& domain)
  : _NPeriodicGhosts(0), _NCellDirectNeib(0),
    _domain(&domain),
    _types(hydro->types),
    _kernrange(hydro->kernp->kernrange)
  { } ;

  NeighbourManager(const ParticleTypeRegister& types, double kernrange,
                    const DomainBox<ndim>& domain)
  : _NPeriodicGhosts(0), _NCellDirectNeib(0),
    _domain(&domain),
    _types(types),
    _kernrange(kernrange)
  { } ;

  void clear() {
    NeighbourManagerDim<ndim>::clear() ;
    neiblist.clear();
    neibdata.clear();
    neib_idx.clear();
    directlist.clear();
  }


  template<class InParticleType>
  void EndSearch(const TreeCellBase<ndim> &cell, const InParticleType* partdata) {
    _EndSearch(cell, partdata, false) ;
  }

  template<class InParticleType>
  void EndSearchGravity(const TreeCellBase<ndim> &cell, const InParticleType* partdata) {
    _EndSearch(cell, partdata, true) ;
  }


  template<class InParticleType>
  void _EndSearch(const TreeCellBase<ndim> &cell, const InParticleType* partdata, bool keep_direct=true) {

    assert(partdata != NULL);


    FLOAT dr[ndim];                      // Relative position vector
    FLOAT drsqd;                         // Distance squared
    FLOAT rc[ndim];                      // Position of cell
    const FLOAT hrangemaxsqd = pow(cell.rmax + _kernrange*cell.hmax,2);
    const FLOAT rmax = cell.rmax;
    for (int k=0; k<ndim; k++) rc[k] = cell.rcell[k];


    const GhostNeighbourFinder<ndim> GhostFinder(*_domain, cell) ;
    int MaxGhosts = GhostFinder.MaxNumGhosts() ;

    assert(!keep_direct || MaxGhosts==1 ||
           (gravcell.size() == 0 && tempdirectneib.size() == 0)) ;

    assert(neibdata.size()==0);
    assert(directlist.size()==0);


    // Now load the particles

    // Start from direct neighbours
    if (keep_direct) {
      for (int ii=0; ii< tempdirectneib.size(); ii++) {
        const int i = tempdirectneib[ii] ;
        const InParticleType& part = partdata[i];
        // Forget immediately: direct particles and particles that do not interact gravitationally
        if (part.flags.is_dead()) continue;
        if (!_types.gravmask[part.ptype]) continue;

        // Now create the particle
        GhostFinder.ConstructGhostsScatterGather(part, neibdata);
        directlist.push_back(neibdata.size()-1);
        neib_idx.push_back(i) ;
      }
    }

    // Now look at the hydro candidate neighbours
    // First the ones that need ghosts to be created on the fly
    int Nneib = directlist.size();
    for (int ii=0; ii < tempperneib.size(); ii++) {
      const int i = tempperneib[ii] ;
      if (partdata[i].flags.is_dead()) continue;

      GhostFinder.ConstructGhostsScatterGather(partdata[i], neibdata);

      while (Nneib < neibdata.size()) {
        int Nmax = neibdata.size() ;
        for (int k=0; k<ndim; k++) dr[k] = neibdata[Nneib].r[k] - rc[k];
        drsqd = DotProduct(dr, dr, ndim);
        FLOAT h2 = rmax + _kernrange*neibdata[Nneib].h;
        if (drsqd < hrangemaxsqd || drsqd < h2*h2) {
          neiblist.push_back(Nneib);
          neib_idx.push_back(i) ;
          Nneib++ ;
        }
        else if (keep_direct) {
          directlist.push_back(Nneib);
          neib_idx.push_back(i) ;
          Nneib++ ;
        }
        else if (Nmax > Nneib) {
          Nmax-- ;
          if (Nmax > Nneib)
            neibdata[Nneib] = neibdata[Nmax];
          neibdata.resize(neibdata.size()-1);
        }
      }// Loop over Ghosts
    }

    // Store the number of periodic particles in the neiblist
    _NPeriodicGhosts = neiblist.size() ;

    // Find those particles that do not need ghosts on the fly
    for (int ii=0; ii<tempneib.size(); ii++) {
      const int i = tempneib[ii];
      if (partdata[i].flags.is_dead()) continue;


      for (int k=0; k<ndim; k++) dr[k] = partdata[i].r[k] - rc[k];
      drsqd = DotProduct(dr,dr,ndim);
      FLOAT h = rmax + _kernrange*partdata[i].h;
      if (drsqd < hrangemaxsqd || drsqd < h*h) {
        neibdata.push_back(partdata[i]);
        neiblist.push_back(Nneib);
        neib_idx.push_back(i) ;
        Nneib++;
      } else if (keep_direct) {
        // Hydro candidates that fail the test get demoted to direct neighbours
        neibdata.push_back(partdata[i]);
        directlist.push_back(Nneib);
        neib_idx.push_back(i) ;
        Nneib++;
      }
    }

    _NCellDirectNeib = directlist.size();
    assert(neibdata.size() == (neiblist.size() + directlist.size()));
  }


  int GetNumAllNeib() {
    return neib_idx.size();
  }

  std::pair<int,ParticleType*> GetNeibI(int i) {
    return make_pair(neib_idx[i],&neibdata[i]);
  }

  template<class InParticleType>
  int GetParticleNeib(const InParticleType& p,const Typemask& hydromask, int** neiblist_p,
      ParticleType** neibdata_p, const bool do_pair_once) {


    FLOAT rp[ndim];                                // Local copy of particle position
    FLOAT draux[ndim];                             // Aux. relative position vector
    FLOAT drsqd;                                   // Distance squared


    for (int k=0; k<ndim; k++) rp[k] = p.r[k];

    const FLOAT hrangesqdi = p.hrangesqd;

    const GhostNeighbourFinder<ndim> GhostFinder(*_domain);

    culled_neiblist.clear();

    // Compute distances and the inverse between the current particle and all neighbours here,
    // for both gather and inactive scatter neibs.
    //-----------------------------------------------------------------------------------------
    for (int j=0; j<neiblist.size(); j++) {
      int jj = neiblist[j] ;

      // Skip if (i) neighbour particle type does not interact hydrodynamically with particle,
      // (ii) neighbour is a dead (e.g. accreted) particle. Additionally, if do_pair_once is true,
      /// (iii) same i.d. as current active particle, (iv) neighbour is on lower timestep level (i.e.
      // timestep is shorter), or (v) neighbour is on same level as current particle but has larger id. value
      // (to only calculate each pair once).
      if (hydromask[neibdata[jj].ptype] == false) continue;

      //if (p.iorig==neibdata[jj].iorig) continue;

      if (do_pair_once) {

        const bool need_interaction =
            (neibdata[jj].flags.is_mirror()) ||
            (p.level > neibdata[jj].level) ||
            (p.level == neibdata[jj].level &&
                p.iorig < neibdata[jj].iorig) ;

        if (not need_interaction) continue;
      }

      for (int k=0; k<ndim; k++) draux[k] = neibdata[jj].r[k] - rp[k];
      if (jj<_NPeriodicGhosts) {
        GhostFinder.NearestPeriodicVector(draux);
      }
      drsqd = DotProduct(draux,draux,ndim)+small_number;

      //          if (drsqd <= small_number) continue ;

      // Compute relative position and distance quantities for pair
      if (drsqd <= hrangesqdi || drsqd <= neibdata[jj].hrangesqd) {
        culled_neiblist.push_back(jj);
        if (jj<_NPeriodicGhosts) for (int k=0; k<ndim; k++) neibdata[jj].r[k] = draux[k]+rp[k];
      }
    }

    *neiblist_p=&culled_neiblist[0];
    *neibdata_p=&neibdata[0];



    return culled_neiblist.size();

  }

  template<class InParticleType>
  ListLength GetParticleNeibGravity(const InParticleType& p,const Typemask& hydromask, const Typemask& gravmask, int** neiblist_p,
      int** directlist_p, int** gravlist_p, ParticleType** neibdata_p, const bool do_grav) {

    FLOAT rp[ndim];
    FLOAT draux[ndim];

    for (int k=0; k<ndim; k++) rp[k] = p.r[k];
    const FLOAT hrangesqdi = p.hrangesqd;

    const GhostNeighbourFinder<ndim> GhostFinder(*_domain);

    culled_neiblist.clear();
    gravlist.clear();
    // Particles that are already in the directlist stay there; we just add the ones that were demoted
    directlist.resize(_NCellDirectNeib);

    // Go through the hydro neighbour candidates and check the distance. The ones that are not real neighbours
    // are demoted to the direct list
    for (int jj=0; jj < neiblist.size(); jj++) {
      int ii = neiblist[jj];

      // Compute relative position and distance quantities for pair
      for (int k=0; k<ndim; k++) draux[k] = neibdata[ii].r[k] - rp[k];
      GhostFinder.NearestPeriodicVector(draux);
      const FLOAT drsqd = DotProduct(draux,draux,ndim) + small_number;

      if (drsqd <= small_number) continue ;

      // Record if neighbour is direct-sum or and SPH neighbour.
      // If SPH neighbour, also record max. timestep level for neighbour
      if (drsqd > hrangesqdi && drsqd >= neibdata[ii].hrangesqd && do_grav) {
        if (gravmask[neibdata[ii].ptype]) directlist.push_back(ii);
      }
      else {
        if (hydromask[neibdata[ii].ptype]){
          culled_neiblist.push_back(ii);
          if (jj<_NPeriodicGhosts) for (int k=0; k<ndim; k++) neibdata[ii].r[k] = draux[k]+rp[k];
        }
        else if (gravmask[neibdata[ii].ptype] && do_grav){
          gravlist.push_back(ii);
        }
      }
    }
    //-----------------------------------------------------------------------------------------
    ListLength listlength;
    listlength.Nhydro=culled_neiblist.size();
    listlength.Ndirect=directlist.size();
    listlength.Ngrav=gravlist.size();

    *neiblist_p = &neiblist[0];
    *directlist_p = &directlist[0];
    *gravlist_p = &gravlist[0];
    *neibdata_p = &neibdata[0];

    return listlength;

  }

};


#endif /* NEIGHBOURMANAGER_H_ */
