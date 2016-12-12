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
  vector<ParticleType > neibpart;

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
    neibpart.clear();
    directlist.clear();
  }


  template<class InParticleType>
  void EndSearch(const TreeCellBase<ndim> &cell, const InParticleType* partdata) {

    assert(partdata != NULL);


    FLOAT dr[ndim];                      // Relative position vector
    FLOAT drsqd;                         // Distance squared
    FLOAT rc[ndim];                      // Position of cell
    const FLOAT hrangemaxsqd = pow(cell.rmax + _kernrange*cell.hmax,2);
    const FLOAT rmax = cell.rmax;
    for (int k=0; k<ndim; k++) rc[k] = cell.rcell[k];


    const GhostNeighbourFinder<ndim> GhostFinder(*_domain, cell) ;
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
        FLOAT h2 = rmax + _kernrange*neibpart[Nneib].h;
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
      FLOAT h = rmax + _kernrange*partdata[i].h;
      if (drsqd < hrangemaxsqd || drsqd < h*h) {
        neibpart.push_back(partdata[i]);
        neiblist.push_back(i);
        Nneib++;
      }
    }

    assert(Nneib==neiblist.size());
    assert(neiblist.size()==neibpart.size());
  }

  template<class InParticleType>
  void EndSearchGravity (const TreeCellBase<ndim> &cell, const InParticleType* partdata) {

    assert(partdata != NULL);


    FLOAT dr[ndim];                      // Relative position vector
    FLOAT drsqd;                         // Distance squared
    FLOAT rc[ndim];                      // Position of cell
    const FLOAT hrangemaxsqd = pow(cell.rmax + _kernrange*cell.hmax,2);
    const FLOAT rmax = cell.rmax;
    for (int k=0; k<ndim; k++) rc[k] = cell.rcell[k];


    const GhostNeighbourFinder<ndim> GhostFinder(*_domain, cell) ;
    int MaxGhosts = GhostFinder.MaxNumGhosts() ;
    assert(MaxGhosts==1);

    assert(neibpart.size()==0);

    // Now load the particles

    // Start from direct neighbours
    for (int ii=0; ii< tempdirectneib.size(); ii++) {
      const int i = tempdirectneib[ii] ;
      const InParticleType& part = partdata[i];
      // Forget immediately: direct particles and particles that do not interact gravitationally
      if (part.flags.is_dead()) continue;
      if (!_types.gravmask[part.ptype]) continue;

      // Now create the particle
      GhostFinder.ConstructGhostsScatterGather(part, neibpart);
      directlist.push_back(neibpart.size()-1);
    }

    // Now look at the hydro candidate neighbours
    for (int ii=0; ii < tempperneib.size(); ii++) {
      const int i = tempperneib[ii] ;
      if (partdata[i].flags.is_dead()) continue;


      GhostFinder.ConstructGhostsScatterGather(partdata[i], neibpart);

      const int index_part = neibpart.size()-1;

      for (int k=0; k<ndim; k++) dr[k] = partdata[i].r[k] - rc[k];
      drsqd = DotProduct(dr, dr, ndim);
      FLOAT h2 = rmax + _kernrange*partdata[i].h;
      if (drsqd < hrangemaxsqd || drsqd < h2*h2) {
        neiblist.push_back(index_part);
      }
      else {
        // Hydro candidates that fail the test get demoted to direct neighbours
        directlist.push_back(index_part);
      }

    }

    _NCellDirectNeib = directlist.size();
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

    const GhostNeighbourFinder<ndim> GhostFinder(*_domain);

    culled_neiblist.clear();

    // Compute distances and the inverse between the current particle and all neighbours here,
    // for both gather and inactive scatter neibs.
    //-----------------------------------------------------------------------------------------
    for (int jj=0; jj<neibpart.size(); jj++) {


      // Skip if (i) neighbour particle type does not interact hydrodynamically with particle,
      // (ii) neighbour is a dead (e.g. accreted) particle. Additionally, if do_pair_once is true,
      /// (iii) same i.d. as current active particle, (iv) neighbour is on lower timestep level (i.e.
      // timestep is shorter), or (v) neighbour is on same level as current particle but has larger id. value
      // (to only calculate each pair once).
      if (hydromask[neibpart[jj].ptype] == false) continue;

      //if (p.iorig==neibpart[jj].iorig) continue;

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

  template<class InParticleType>
  ListLength GetParticleNeibGravity(const InParticleType& p,const Typemask& hydromask, const Typemask& gravmask, int** neiblist_p,
      int** directlist_p, int** gravlist_p, ParticleType** neibpart_p, const bool do_grav) {

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
      for (int k=0; k<ndim; k++) draux[k] = neibpart[ii].r[k] - rp[k];
      GhostFinder.NearestPeriodicVector(draux);
      const FLOAT drsqd = DotProduct(draux,draux,ndim) + small_number;

      if (drsqd <= small_number) continue ;

      // Record if neighbour is direct-sum or and SPH neighbour.
      // If SPH neighbour, also record max. timestep level for neighbour
      if (drsqd > hrangesqdi && drsqd >= neibpart[ii].hrangesqd && do_grav) {
        if (gravmask[neibpart[ii].ptype]) directlist.push_back(ii);
      }
      else {
        if (hydromask[neibpart[ii].ptype]){
          culled_neiblist.push_back(ii);
        }
        else if (gravmask[neibpart[ii].ptype] && do_grav){
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
    *neibpart_p = &neibpart[0];

    return listlength;

  }

};


#endif /* NEIGHBOURMANAGER_H_ */
