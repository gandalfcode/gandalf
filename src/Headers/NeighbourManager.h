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

//=================================================================================================
//  Class NeighbourManager
/// \brief   Template class for neighbour searches.
/// \details Handles issues related to positioning of ghosts, along with culling the neighbour list
///          down to only those needed in a tree walk.
/// \author  G. Rosotti, R. A. Booth
/// \date    15/12/2016
//=================================================================================================
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
  vector<int> smoothgravlist;

  vector<FLOAT> dr;
  vector<FLOAT> drmag;

  const DomainBox<ndim>* _domain;
  const ParticleTypeRegister* _types ;
  double _kernrange ;

  typedef struct {}        _true_type ;
  typedef struct {char c;} _false_type ;

public:
  using NeighbourManagerDim<ndim>::tempneib;
  using NeighbourManagerDim<ndim>::tempperneib;
  using NeighbourManagerDim<ndim>::tempdirectneib;
  using NeighbourManagerDim<ndim>::gravcell;

  NeighbourManager(const Hydrodynamics<ndim>* hydro, const DomainBox<ndim>& domain)
  : _NPeriodicGhosts(0), _NCellDirectNeib(0),
    _domain(&domain),
    _types(&(hydro->types)),
    _kernrange(hydro->kernp->kernrange)
  { } ;

  NeighbourManager(const ParticleTypeRegister& types, double kernrange,
      const DomainBox<ndim>& domain)
  : _NPeriodicGhosts(0), _NCellDirectNeib(0),
    _domain(&domain),
    _types(&types),
    _kernrange(kernrange)
  { } ;

  void clear() {
    NeighbourManagerDim<ndim>::clear() ;
    neiblist.clear();
    neibdata.clear();
    neib_idx.clear();
    directlist.clear();
  }

  //===============================================================================================
  // EndSearch
  /// \brief Colelct particle data needed for the neighbours and cull distant particles that do
  ///        not interact with the cell hydrodyanimcally
  //===============================================================================================
  template<class InParticleType>
  void EndSearch(const TreeCellBase<ndim> &cell, const InParticleType* partdata) {
    _EndSearch(cell, partdata, false) ;
  }
  //===============================================================================================
  // EndSearchGravity
  /// \brief Collect particle data needed for the neighbours, and demote particles that do not
  ///        interact hydrodynamically to the directlist.
  //===============================================================================================
  template<class InParticleType>
  void EndSearchGravity(const TreeCellBase<ndim> &cell, const InParticleType* partdata) {
    _EndSearch(cell, partdata, true) ;
  }

  //===============================================================================================
  //  GetNumAllNeib
  /// \brief Return the total number of all particles found in the walk.
  //===============================================================================================
  int GetNumAllNeib() {
    return neib_idx.size();
  }
  //===============================================================================================
  //  GetNeibI
  /// \brief Get the index of thew i-th neighbour in the original particle array, and a pointer to
  ///        to the reduced neighbour data stored.
  //===============================================================================================
  std::pair<int,ParticleType*> GetNeibI(int i) {
    return make_pair(neib_idx[i],&neibdata[i]);
  }

  //===============================================================================================
  //  GetParticleNeib
  /// \brief    Get the list of particles that interact hydrodynamically with the Particle, p.
  /// \details  Gets a trimmed list of particles that interact with hydrodynamically with Particle,
  ///           p. This is determined by whether or not the particles smoothing spheres overlap,
  ///           and hydromask, which specifies the types needed (hydromask[ptype] == true for the
  ///           required particles). Also, if do_pair_once is included then the neighbour is only
  ///           included if this interaction will not have been already calculated in the update
  ///           for the neighbour itself.
  /// \returns  Returns the number of neighbours found. Additionally, neibdata_p is set to the
  ///           list of particles stored by the data manager, and neiblist_p is set to the indices
  ///           within neibdata_p of the required neighbours.
  //===============================================================================================
  template<class InParticleType>
  int GetParticleNeib
  (const InParticleType& p,                            ///< [in] Particle to collect the neibs for
   const Typemask& hydromask,                          ///< [in] Boolean flags listing types we need
   int** neiblist_p,                                   ///< [out] List of particles needed.
   ParticleType** neibdata_p,                          ///< [out] List of particle data
   const bool do_pair_once)                            ///< [in]
  {
    if (do_pair_once)
      TrimNeighbourLists<InParticleType,_true_type>(p, hydromask, false) ;
    else
      TrimNeighbourLists<InParticleType,_false_type>(p, hydromask, false) ;


    *neiblist_p=&culled_neiblist[0];
    *neibdata_p=&neibdata[0];
    return culled_neiblist.size();

  }

  //===============================================================================================
  //  GetParticleNeibGravity
  /// \brief    Get the list of particles that interact with the Particle, p.
  /// \details  As with GetParticleNeib, this function returns the list of neighbours that interact
  ///           hydrodynamically with the Particle, p. Additionally this function also generates
  ///           the list of particles needed for gravitational interactions, including the list of
  ///           smoothed and unsmoothed contributions.
  /// \returns  Returns the number of neighbours found. Additionally, neibdata_p is set to the
  ///           list of particles stored by the data manager, and neiblist_p, directlist_p and
  ///           smoothgravlist_p are set to the indices within neibdata_p of the required
  ///           neighbours.
  //===============================================================================================
  template<class InParticleType>
  ListLength GetParticleNeibGravity
  (const InParticleType& p,                            ///< [in] Particle to collect the neibs for
   const Typemask& hydromask,                          ///< [in] Type flags for hydro interactions
   int** neiblist_p,                                   ///< [out] List of hydro neigbour indices
   int** directlist_p,                                 ///< [out] List of unsmoothed gravity neibs
   int** smoothgravlist_p,                             ///< [out] List of smoothed gravity neibs
   ParticleType** neibdata_p) {                        ///< [out] List of particle data.

    TrimNeighbourLists<InParticleType,_false_type>(p, hydromask, true) ;

    ListLength listlength;
    listlength.Nhydro=culled_neiblist.size();
    listlength.Ndirect=directlist.size();
    listlength.Ngrav=smoothgravlist.size();

    *neiblist_p = &neiblist[0];
    *directlist_p = &directlist[0];
    *smoothgravlist_p = &smoothgravlist[0];
    *neibdata_p = &neibdata[0];

    return listlength;

  }

private:

  //===============================================================================================
  //  _EndSearch
  /// \details Gathers the particle data for the cell, and trims the list of hydrodynamic particles.
  ///          If keep_direct is true, then these become direct-list gravitational neighbours,
  ///          otherwise they are discarded. Periodic corrections are included.
  //===============================================================================================
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

     Typemask gravmask = _types->gravmask ;

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
         if (!gravmask[part.ptype]) continue;

         // Now create the particle
         GhostFinder.ConstructGhostsScatterGather(part, neibdata);
         directlist.push_back(neibdata.size()-1);
         neib_idx.push_back(i) ;
       }
     }

     assert(directlist.size() == neibdata.size() &&
         neib_idx.size()   == neibdata.size()) ;

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
         else if (keep_direct && gravmask[neibdata[Nneib].ptype]) {
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
       } else if (keep_direct && gravmask[neibdata[Nneib].ptype]) {
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


  //===============================================================================================
  //  TrimNeighbourLists
  /// \detail This function trims the neighbour lists for a given particle, determining whether
  ///         neighbours are needed for hydro or gravity. Periodic corrections are applied.
  //===============================================================================================
  template<class InParticleType, class do_pair_once>
  void TrimNeighbourLists(const InParticleType& p, const Typemask& hydromask, bool keep_grav)
  {
    FLOAT rp[ndim];
    FLOAT draux[ndim];

    for (int k=0; k<ndim; k++) rp[k] = p.r[k];
    const FLOAT hrangesqdi = p.hrangesqd;

    const GhostNeighbourFinder<ndim> GhostFinder(*_domain);

    Typemask gravmask = _types->gravmask ;

    culled_neiblist.clear();
    smoothgravlist.clear();
    // Particles that are already in the directlist stay there; we just add the ones that were demoted
    directlist.resize(_NCellDirectNeib);

    // Go through the hydro neighbour candidates and check the distance. The ones that are not real neighbours
    // are demoted to the direct list
    for (int jj=0; jj < neiblist.size(); jj++) {
      int ii = neiblist[jj];

      // If do_pair_once is true then only get the neighbour for the first of the two times the
      if (not _first_appearance(p, neibdata[ii], do_pair_once())) continue ;

      // Compute relative position and distance quantities for pair
      for (int k=0; k<ndim; k++) draux[k] = neibdata[ii].r[k] - rp[k];
      if (jj < _NPeriodicGhosts && hydromask[neibdata[ii].ptype])
        GhostFinder.ApplyPeriodicDistanceCorrection(neibdata[ii].r, draux);

      const FLOAT drsqd = DotProduct(draux,draux,ndim) + small_number;

      //if (drsqd <= small_number) continue ;

      // Record if neighbour is direct-sum or and SPH neighbour.
      // If SPH neighbour, also record max. timestep level for neighbour
      if (drsqd > hrangesqdi && drsqd >= neibdata[ii].hrangesqd) {
        if (keep_grav && gravmask[neibdata[ii].ptype]) directlist.push_back(ii);
      }
      else {
        if (hydromask[neibdata[ii].ptype]){
          culled_neiblist.push_back(ii);
        }
        else if (keep_grav && gravmask[neibdata[ii].ptype]) {
          smoothgravlist.push_back(ii);
        }
      }
    }
  }

  template<class InParticleType>
  bool _first_appearance(const InParticleType& p, const ParticleType& neibpart, _true_type) {
    // Get the neighbour only if this is it's first appearance, which will be true if any of
    // the following criteria are met.
    // (i)   the particle is a ghost
    // (ii)  same i.d. as current active particle
    // (iii) neighbour is on lower timestep level
    // (iv)  the neighbour is on same timestep level as current particle but has larger id value
    return
        (neibpart.flags.is_mirror()) ||
        (p.level > neibpart.level) ||
        (p.level == neibpart.level &&  p.iorig < neibpart.iorig) ;
  }

  template<class InParticleType>
  bool _first_appearance(const InParticleType& p, const ParticleType& neibpart, _false_type) {
    return true ;
  }


};


#endif /* NEIGHBOURMANAGER_H_ */
