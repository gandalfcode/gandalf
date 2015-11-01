//=================================================================================================
//  Dust.cpp
//  Contains the main implementation for the dust drag force classes.
//  This file is part of GANDALF :
//  Graphical Astrophysics code for N-body Dynamics And Lagrangian Fluids
//  https://github.com/gandalfcode/gandalf
//  Contact : gandalfcode@gmail.com
//
//  Copyright (C) 2013  D. A. Hubber, G. Rosotti
//
//  GANDALF is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 2 of the License, or
//  (at your option) any later version.
//
//  GANDALF is distributed in the hope that it will be useful, but
//  WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  General Public License (http://www.gnu.org/licenses) for more details.
//=================================================================================================

#include <assert.h>
#include <cmath>
#include <vector>
#include "CodeTiming.h"
#include "Debug.h"
#include "DragLaws.h"
#include "Dust.h"
#include "Particle.h"
#include "Precision.h"
#include "Tree.h"

//=================================================================================================
//  Class DustSphNgbFinder
/// \brief   DustSphNgbFinder class definition.
/// \details Base class for drag forces, finds neighbours for drag force. Implementer must provide
///          the routines for computing the drag force.
/// \author  R. A. Booth
/// \date    17/10/2015
//=================================================================================================
template<int ndim, template<int> class ParticleType >
class DustSphNgbFinder
: public DustBase<ndim>
{
public:
	DustSphNgbFinder(TreeBase<ndim> * t,TreeBase<ndim> * gt=NULL)
	: _tree(t), _ghosttree(gt)
	{ } ;
	virtual ~DustSphNgbFinder() { } ;

#ifdef MPI_PARALLEL
	void set_mpi_tree(TreeBase<ndim>* t)
	{ mpighosttree = t ; }
#endif

protected:
	template<class Interp>
	void FindNeibAndDoInterp(int ,int, ParticleType<ndim> *,Typemask, Interp&) ;

	template<class InterpData>
	void FindNeibForces() ;
private:

	TreeBase<ndim>* _tree, *_ghosttree ;   ///< Pointer to neighbour tree
#if defined MPI_PARALLEL
	TreeBase<ndim>* mpighosttree;        ///< Pointer to pruned tree arrays
#endif
protected:

};

/*
//=================================================================================================
//  Class DustFull
/// \brief   DustFull class definition.
/// \details Computes the drag force on the dust and it's back-reaction on the gas
/// \author  R. A. Booth
/// \date    17/10/2015
//=================================================================================================
template<int ndim, template<int> class ParticleType, class StoppingTime, class Kernel>
class DustFull
: public DustSphNgbFinder<ndim,ParticleType<ndim> >
{
public:
	DustFull() {} ;
	void UpdateAllDragForces(int NPart, Ntot, Particle<ndim> *sph_gen) ;
private:
	void ComputeDragForces(const int i,
						   const int nNeib, const int *neiblist,
						   const FLOAT *drmag, const FLOAT *dr,
				           Particle<ndim> &part, Particle<ndim> *neibpart_gen) ;
};
*/

//=================================================================================================
//  Class DustTestPartInterp
/// \brief   DustTestPartInterp class definition.
/// \details Collector class that grabs the required interpolation data for test particle drag
///          force calculation.
/// \author  R. A. Booth
/// \date    17/10/2015
//=================================================================================================
template<int ndim, template <int> class ParticleType, class StoppingTime, class Kernel>
class DustInterpolant
{
  struct DustTestPartInterp
  {
	DustTestPartInterp() { } ;
	DustTestPartInterp(const ParticleType<ndim>& p){
	  cs = p.sound ;
	  for (int k=0; k < ndim; ++k){
		  v[k] = p.v[k] ;
		  a[k] = p.a[k] ;
	  }
	}

	FLOAT cs;
	FLOAT v[ndim] ;
	FLOAT a[ndim] ;
  };
public:

  DustInterpolant(const StoppingTime& ts,const Kernel &k, FLOAT h_factor, FLOAT h_conv )
   : kern(k), t_stop(ts), h_fac(h_factor), h_converge(h_conv)
  { } ;

  typedef DustTestPartInterp DataType;

  int DoInterpolate(const int, const int , FLOAT hmax,
                    const std::vector<FLOAT>&,const std::vector<FLOAT>&,
                    const std::vector<DataType>&, ParticleType<ndim>&) ;

private:
  Kernel kern ;
  StoppingTime t_stop ;

  FLOAT h_fac ;
  FLOAT h_converge ;
};
//=================================================================================================
//  Class DustTestParticle
/// \brief   DustTestParticle class definition.
/// \details Compute the drag forces on the dust in the test particle limit.
/// \author  R. A. Booth
/// \date    17/10/2015
//=================================================================================================
template<int ndim, template <int> class ParticleType, class StoppingTime, class Kernel>
class DustTestParticle
: public DustSphNgbFinder<ndim,ParticleType>
{
	using DustSphNgbFinder<ndim,ParticleType>::FindNeibAndDoInterp ;
public:

	typedef DustInterpolant<ndim, ParticleType, StoppingTime, Kernel> DI ;

	DustTestParticle(DI Interp, TreeBase<ndim> * t, TreeBase<ndim> * gt=NULL)
	 :  DustSphNgbFinder<ndim, ParticleType>(t, gt),
	    _interp(Interp)
    { } ;
	void UpdateAllDragForces(int NPart, int Ntot, Particle<ndim> *sph_gen)
	{
		debug2("[DustTestParticle::UpdateAllDragForces]") ;

		ParticleType<ndim>* sphdata = reinterpret_cast<ParticleType<ndim>*>(sph_gen) ;

		Typemask mask = {false} ;
		mask[gas_type] = true ;

		FindNeibAndDoInterp(NPart, Ntot, sphdata, mask, _interp) ;

		//TODO: Check: Should this be NPart or Ntot?
		for (int i(0); i < NPart; ++i)
			if (sphdata[i].active && sphdata[i].ptype == dust_type){
			    for (int k(0); k < ndim; ++k)
			        sphdata[i].a[k] += sphdata[i].a_dust[k] ;
			}
	}
private:
  DI _interp ;
};


//=================================================================================================
//  Member function DustSphNgbFinder::FindNgbINterpolate
// \brief    Driver routine for interpolation to location of dust particles
/// \details This function finds the neighbours for a given dust particle that are within
///          a smoothing length to be used for interpolating gas quantities to the location of dust
///          particles. Simultaneously finds the smoothing length.
/// \author  R. A. Booth
/// \date    20/10/2015
//=================================================================================================
template<int ndim, template<int> class ParticleType>
template<class Interpolant>
void DustSphNgbFinder<ndim, ParticleType>::FindNeibAndDoInterp
(int Nhydro,                               ///< [in] No. of SPH particles
 int Ntot,                                 ///< [in] No. of SPH + ghost particles
 ParticleType<ndim> *sphdata,              ///< [inout] Pointer to SPH ptcl array
 Typemask mask,                            ///< [in] Types to include in the interpolation
 Interpolant& Interp)                      ///< [in] Interpolating function
{
  using std::vector ;

  typedef typename Interpolant::DataType InterpData;


  int cactive;                             // No. of active tree cells
  TreeCellBase<ndim> **celllist;           // List of active tree cells

#ifdef MPI_PARALLEL
  int Nactivetot = 0;                      // Total number of active particles
  //double twork = timing->WallClockTime();  // Start time (for load balancing)
#endif

  debug2("[DustSphNgbFinder::FindNeibAndDoInterp]");
  //timing->StartTimingSection("DUST_GAS_INTERPOLATE");

  // Find list of all cells that contain active particles
  celllist = new TreeCellBase<ndim>*[_tree->MaxNumCells()];
  cactive = _tree->ComputeActiveCellPointers(celllist);
  assert(cactive <= _tree->MaxNumCells());


  // Set-up all OMP threads
  //===============================================================================================
#pragma omp parallel default(none) shared(cactive,celllist,cout,nbody,sphdata) private()
  {
#if defined _OPENMP
    const int ithread = omp_get_thread_num();
#else
    const int ithread = 0;
#endif
    int celldone;                              // Flag if cell is done
    int cc;                                    // Aux. cell counter
    int i;                                     // Particle id
    int j;                                     // Aux. particle counter
    int jj;                                    // Aux. particle counter
    int k;                                     // Dimension counter
    int Nactive;                               // No. of active particles in cell
    int Ngather;                               // No. of gather neighbours
    int Nneib;                                 // No. of neighbours from tree-walk
    int okflag;                                // Flag if particle is done
    FLOAT draux[ndim];                         // Aux. relative position vector var
    FLOAT drsqdaux;                            // Distance squared
    FLOAT hrangesqd;                           // Kernel extent
    FLOAT hmax;                                // Maximum smoothing length
    FLOAT rp[ndim];                            // Local copy of particle position
    int Nneibmax = 2000 ;

    vector<int>        neiblist(Nneibmax);     // Local array of neighbour particle ids
    vector<FLOAT>      pos(ndim*Nneibmax);     // Local reduced array of neighbour potentials
    vector<FLOAT>      drsqd(Nneibmax);        // Local array of distances (squared)
    vector<FLOAT>      m(Nneibmax);            // Local array of particle masses
    vector<FLOAT>      m2(Nneibmax);           // Local array of particle masses (reduced)
    vector<int>        ptype(Nneibmax);        // Local array of particle types
    vector<InterpData> data(Nneibmax);         // Local array of data to be interpolated
    vector<InterpData> data2(Nneibmax);        // Local array of data to be interpolated (reduced)

    vector<int> activelist(_tree->MaxNumPartInLeafCell()) ;  // Local array of active particles ids
    vector<ParticleType<ndim> > activepart(_tree->MaxNumPartInLeafCell()) ;  // Local array of active particles

    FLOAT kernrangesqd = pow(_tree->MaxKernelRange(),2) ;


    // Loop over all active cells
    //=============================================================================================
#pragma omp for schedule(guided)
    for (cc=0; cc<cactive; cc++) {
      TreeCellBase<ndim>& cell = *(celllist[cc]);
      celldone = 1;
      hmax = cell.hmax;


      // If hmax is too small so the neighbour lists are invalid, make hmax
      // larger and then recompute for the current active cell.
      //-------------------------------------------------------------------------------------------
      do {
        hmax = (FLOAT) 1.05*hmax;
        celldone = 1;

        // Find list of active particles in current cell
        Nactive = _tree->ComputeActiveParticleList(cell, sphdata, &(activelist[0]));
        for (j=0; j<Nactive; j++) {
        	activepart[j] = sphdata[activelist[j]];

        	// Make sure that the smoothing length is at least as big as the dust one
        	if (activepart[j].ptype == dust_type && 1.05 * activepart[j].h_dust > hmax)
        		hmax = 1.05 * activepart[j].h_dust ;
        }

        // Compute neighbour list for cell from particles on all trees
        do {
          // TODO: These functions should accept vectors so that the functions can
          //       automatically manage their size.
          Nneib = 0;
          Nneib = _tree->ComputeGatherNeighbourList(cell,sphdata,hmax,Nneibmax,Nneib,&(neiblist[0]));
          Nneib = _ghosttree->ComputeGatherNeighbourList(cell,sphdata,hmax,Nneibmax,Nneib,
        												 &(neiblist[0]));
#ifdef MPI_PARALLEL
          Nneib = mpighosttree->ComputeGatherNeighbourList(cell,sphdata,hmax,Nneibmax,Nneib,
        				                                   &(neiblist[0]));
#endif


          // If there are too many neighbours so the buffers are filled,
          // reallocate the arrays and recompute the neighbour lists.
          if (Nneib < 0) {
            Nneibmax *= 2;

            neiblist.resize(Nneibmax);;
            drsqd.resize(Nneibmax);
            m.resize(Nneibmax);
            m2.resize(Nneibmax);
            pos.resize(Nneibmax*ndim);
            ptype.resize(Nneibmax);
            data.resize(Nneibmax);
            data2.resize(Nneibmax);
          }
        }  while (Nneib < 0) ;

        // Make local copies of important neib information (mass and position)
        for (jj=0; jj<Nneib; jj++) {
          j        = neiblist[jj];

          m[jj]     = sphdata[j].m;
          ptype[jj] = sphdata[j].ptype;
          for (k=0; k<ndim; k++) pos[ndim*jj + k] = sphdata[j].r[k];

          data[jj] = InterpData(sphdata[j]) ;
        }

        // Loop over all active particles in the cell
        //-----------------------------------------------------------------------------------------
        for (j=0; j<Nactive; j++) {
          i = activelist[j];

          // Skip non-dust particles
          if (activepart[j].ptype != dust_type)
        	 continue ;

          for (k=0; k<ndim; k++) rp[k] = activepart[j].r[k];

          // Set gather range as current h multiplied by some tolerance factor
          hrangesqd = kernrangesqd*hmax*hmax;
          Ngather = 0;

          // Compute distance (squared) to all
          //---------------------------------------------------------------------------------------
          for (jj=0; jj<Nneib; jj++) {
        	// Only include particles of appropriate types in density calculation
        	if (!mask[ptype[jj]]) continue ;

            for (k=0; k<ndim; k++) draux[k] = pos[ndim*jj + k] - rp[k];
            drsqdaux = DotProduct(draux,draux,ndim) + small_number;

            // Record distance squared and masses for all potential gather neighbours
            if (drsqdaux <= hrangesqd) {
              drsqd[Ngather] = drsqdaux;
              m2[Ngather]    = m[jj];
              data2[Ngather] = data[jj] ;
              Ngather++;
            }

          }
          //---------------------------------------------------------------------------------------

          // Validate that gather neighbour list is correct
#if defined(VERIFY_ALL)
          if (neibcheck) this->CheckValidNeighbourList(i, Ntot, Nneib, neiblist, sphdata, "gather");
#endif

          // Compute smoothing length and other gather properties for ptcl i
          okflag = Interp.DoInterpolate(i, Ngather, hmax, m2, drsqd, data2, activepart[j]) ;

      	// If h-computation is invalid, then break from loop and recompute
          // larger neighbour lists
          if (okflag == 0) {
            celldone = 0;
            break;
          }

        }
        //-----------------------------------------------------------------------------------------

      } while (celldone == 0);
      //-------------------------------------------------------------------------------------------

      // Once cell is finished, copy all active particles back to main memory
      for (j=0; j<Nactive; j++) sphdata[activelist[j]] = activepart[j];
    }
  }
  //===============================================================================================

  // Compute time spent in routine and in each cell for load balancing
  /*
#ifdef MPI_PARALLEL
  twork = timing->WallClockTime() - twork;
  for (int cc=0; cc<cactive; cc++) Nactivetot += celllist[cc].Nactive;
  for (int cc=0; cc<cactive; cc++) {
    int c = celllist[cc].id;
    tree->celldata[c].worktot += twork*(DOUBLE) tree->celldata[c].Nactive / (DOUBLE) Nactivetot;
  }
  cout << "Time computing dust smoothing lengths : " << twork << "     Nactivetot : " << Nactivetot << endl;
#endif
*/

  //timing->EndTimingSection("DUST_GAS_INTERPOLATE");

  return;
}

/*
//=================================================================================================
//  Member function DustSphNgbFinder::FindNgbForces
// \brief    Driver routine for pairwise particle forces
/// \details This function finds the pairs of particles that are involved in a given force
///          calculation.
/// \author  R. A. Booth
/// \date    20/10/2015
//=================================================================================================
template<int ndim, template<int> class ParticleType, template<int> class TreeCell>
template<class InterpData, class InterpFuction>
void DustSphNgbFinder::FindNeibForces
(int Nhydro,                              ///< [in] No. of SPH particles
  int Ntot,                                ///< [in] No. of SPH + ghost particles
  Particle<ndim> *sph_gen,              ///< [inout] Pointer to SPH ptcl array
  Sph<ndim> *sph,                          ///< [in] Pointer to SPH object
  Nbody<ndim> *nbody)                      ///< [in] Pointer to N-body object
{
  using std::vector ;

  int cactive;                             // No. of active cells
  TreeCell<ndim> *celllist;                // List of active tree cells
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);
#ifdef MPI_PARALLEL
  int Nactivetot = 0;                      // Total number of active particles
  double twork = timing->WallClockTime();  // Start time (for load balancing)
#endif

  debug2("[GradhSphTree::UpdateAllSphHydroForces]");
  timing->StartTimingSection("DUST_GAS_FORCES");


  // Find list of all cells that contain active particles
#if defined (MPI_PARALLEL)
  celllist = new TreeCell<ndim>[tree->Ncellmax];
#else
  celllist = new TreeCell<ndim>[tree->];
#endif
  cactive = tree->ComputeActiveCellList(celllist);

  // If there are no active cells, return to main loop
  if (cactive == 0) {
    timing->EndTimingSection("DUST_GAS_FORCES");
    return;
  }

  // Update ghost tree smoothing length values here
  tree->UpdateAllHmaxValues(sphdata);
  if (ghosttree->Ntot > 0) ghosttree->UpdateAllHmaxValues(sphdata);


  // Set-up all OMP threads
  //===============================================================================================
#pragma omp parallel default(none) shared(cactive,celllist,nbody,sph,sphdata)
  {
#if defined _OPENMP
    const int ithread = omp_get_thread_num();
#else
    const int ithread = 0;
#endif
    int cc;                                      // Aux. cell counter
    int i;                                       // Particle id
    int j;                                       // Aux. particle counter
    int jj;                                      // Aux. particle counter
    int k;                                       // Dimension counter
    int Nactive;                                 // ..
    int Nneib;                                   // ..
    int Nhydroaux;                               // ..
    FLOAT draux[ndim];                           // Aux. relative position vector
    FLOAT drsqd;                                 // Distance squared
    FLOAT hrangesqdi;                            // Kernel gather extent
    FLOAT rp[ndim];                              // Local copy of particle position
    Typemask hydromask;                          // Mask for computing hydro forces
    int Nneibmax    = Nneibmaxbuf[ithread];      // ..
    vector<int> * activelist = activelistbuf[ithread];    // ..
    vector<int>       neiblist(Nneibmax);        // Local array of neighbour particle ids
    vector<int>       sphlist(Nneibmax) ;       // ..
    vector<FLOAT>     dr(ndim*Nneibmax);         // Local reduced array of neighbour potentials
    vector<FLOAT>     rmag(Nneibmax);            // Local array of distances (squared)
    vector<FLOAT>     invdrmag(Nneibmax);        // ...
    vector<ParticleType<ndim>> activepart(Nneibmax);  // Local array of particles
    vector<ParticleType<ndim>> neibpart(Nneibmax);    // Local array of particles (reduced)

    // Loop over all active cells
    //=============================================================================================
#pragma omp for schedule(guided)
    for (cc=0; cc<cactive; cc++) {
      TreeCell<ndim>& cell = celllist[cc];

      // Find list of active particles in current cell
      Nactive = tree->ComputeActiveParticleList(cell,sphdata, &(activelist[0]));

      // Make local copies of active particles
      for (j=0; j<Nactive; j++) {
        activepart[j] = sphdata[activelist[j]];
      }


      // Compute neighbour list for cell from real and periodic ghost particles
      do {
    	// TODO: These functions should accept vectors so that the functions can
    	//       automatically manage their size.
        Nneib = 0;
        Nneib = tree->ComputeNeighbourList(cell, sphdata, Nneibmax, Nneib,
        								   &(neiblist[0]), &(neibpart[0]));
        Nneib = ghosttree->ComputeNeighbourList(cell, sphdata, Nneibmax, Nneib,
        										&(neiblist[0]), &(neibpart[0]));

        // If there are too many neighbours, reallocate the arrays and
        // recompute the neighbour list.

        if (Nneib == -1) {
    	  Nneibmax *= 2 ;

    	  neibpart.resize(Nneibmax);
    	  invdrmag.resize(Nneibmax);
    	  drmag.resize(Nneibmax);
    	  dr.resize(ndim*Nneibmax);
    	  sphlist.resize(Nneibmax);
    	  neiblist.resize(Nneibmax);
        }
      }
      while (Nneib == -1) ;

      // Loop over all active particles in the cell
      //-------------------------------------------------------------------------------------------
      for (j=0; j<Nactive; j++) {
        i = activelist[j];

        for (k=0; k<ndim; k++) rp[k] = activepart[j].r[k];
        for (k=0; k<Nhydrotypes; k++) dragmask[k] = sph->types[activepart[j].ptype].hydromask[k];
        hrangesqdi = activepart[j].hrangesqd;
        Nhydroaux = 0;

        // Validate that gather neighbour list is correct
#if defined(VERIFY_ALL)
        if (neibcheck) this->CheckValidNeighbourList(i, sph->Nhydro + sph->NPeriodicGhost,
                                                     Nneib, neiblist, sphdata, "all");
#endif

        // Compute distances and the inverse between the current particle and all neighbours here,
        // for both gather and inactive scatter neibs.  Only consider particles with j > i to
        // compute pair forces once unless particle j is inactive.
        //-----------------------------------------------------------------------------------------
        for (jj=0; jj<Nneib; jj++) {

          // Skip current active particle or if neighbour type is not active for hydro forces
          if (neiblist[jj] == i || dragmask[neibpart[jj].ptype] == false) continue;

          for (k=0; k<ndim; k++) draux[k] = neibpart[jj].r[k] - rp[k];
          drsqd = DotProduct(draux, draux, ndim) + small_number;

          // Compute relative position and distance quantities for pair
          if (drsqd <= hrangesqdi || drsqd <= neibpart[jj].hrangesqd) {
            drmag[Nhydroaux] = sqrt(drsqd);
            invdrmag[Nhydroaux] = (FLOAT) 1.0/drmag[Nhydroaux];
            for (k=0; k<ndim; k++) dr[Nhydroaux*ndim + k] = draux[k]*invdrmag[Nhydroaux];
            sphlist[Nhydroaux] = jj;
            Nhydroaux++;
          }

        }
        //-----------------------------------------------------------------------------------------

        // Compute all neighbour contributions to the drag forces
        sph->ComputeSphHydroForces(i,Nhydroaux,sphlist,drmag,invdrmag,dr,activepart[j],neibpart);

      }
      //-------------------------------------------------------------------------------------------


      // Add all active particles contributions to main array
      for (j=0; j<Nactive; j++) {
        i = activelist[j];
        for (k=0; k<ndim; k++) sphdata[i].a[k] += activepart[j].a_dust[k] ;
        sphdata[i].dudt = activepart[j].dudt ;
      }

    }
  }
  //===============================================================================================

  // Compute time spent in routine and in each cell for load balancing
#ifdef MPI_PARALLEL
  twork = timing->WallClockTime() - twork;
  for (int cc=0; cc<cactive; cc++) Nactivetot += celllist[cc].Nactive;
  for (int cc=0; cc<cactive; cc++) {
    int c = celllist[cc].id;
    tree->celldata[c].worktot += twork*(DOUBLE) tree->celldata[c].Nactive / (DOUBLE) Nactivetot;
  }
  cout << "Time computing dust forces : " << twork << "     Nactivetot : " << Nactivetot << endl;
#endif


  timing->EndTimingSection("DUST_GAS_FORCES");

  return;
}



template<int ndim, class StoppingTime, class Kernel>
void DustFull<ndim, StoppingTime, Kernel>::ComputeDragForces
(const int i,                         ///< [in] id of particle
 const int Nneib,                     ///< [in] No. of neins in neibpart array
 const int *neiblist,                 ///< [in] id of gather neibs in neibpart
 const FLOAT *drmag,                  ///< [in] Distances of gather neighbours
 const FLOAT *dr,                     ///< [in] Position vector of gather neibs
 Particle<ndim> &part,                ///< [inout] Particle i data
 Particle<ndim> *neibpart_gen)        ///< [inout] Neighbour particle data
{
  int j;                               // Neighbour list id
  int jj;                              // Aux. neighbour counter
  int k;                               // Dimension counter
  FLOAT draux[ndim];                   // Relative position vector
  FLOAT dv[ndim];                      // Relative velocity vector
  FLOAT da[ndim];                      // Relative acceleration vector
  FLOAT dvdr;                          // Dot product of dv and dr
  FLOAT dadr;                          // Dot product of da and dr
  FLOAT wkern;                         // Value of w1 kernel function for part i
  FLOAT d2g ;                          // Dust to gas ratio
  FLOAT t_s ;                          // Stopping time
  FLOAT Xi ;                           // Stopping factor
  FLOAT S ;                            // Drag term
  GradhSphParticle<ndim>& parti = static_cast<GradhSphParticle<ndim>& > (part);
  GradhSphParticle<ndim>* neibpart = static_cast<GradhSphParticle<ndim>* > (neibpart_gen);

  // Some basic sanity-checking in case of invalid input into routine
  assert(neibpart_gen != NULL);
  assert(parti.itype != dead);

  // Loop over all potential neighbours in the list
  //-----------------------------------------------------------------------------------------------
  for (jj=0; jj<Nneib; jj++) {
    j = neiblist[jj];
    assert(neibpart[j].itype != dead);

    if (part[i].ptype == gas_type)
    	wkern = parti.hfactor*kern.wdrag(drmag[jj]*parti.invh);
    else
    	wkern = neibpart[j].hfactor*kern.wdrag(drmag[jj]*neibpart[j].invh);

    wkern *= neibpart[j].m / neibpart[j].rho ;

    for (k=0; k<ndim; k++)
    {
    	draux[k] = dr[jj*ndim + k];
    	dv[k] = parti.v[k] - neibpart[j].v[k] ;
    	da[k] = parti.a[k] - neibpart[j].a[k] ;
    }
    dvdr = DotProduct(draux, dv, ndim);
    dadr = DotProduct(draux, da, ndim);


    //---------------------------------------------------------------------------------------------
    // Compute stopping time
    if (parti.ptype == gas_type)
    {
    	d2g = neibpart[j].rho/part[i].rho ;
    	t_s = StoppingTime(parti, d2g) ;
    }
    else
    {
    	d2g = part[i].rho/neibpart[j].rho ;
    	t_s = StoppingTime(neibpart[j], d2g) ;
    }

    //---------------------------------------------------------------------------------------------
    // Evaluate the drag term
    Xi = (1 - exp(- parti.dt/t_s)) / (1 + d2g) ;
    S = dvdr * Xi - dadr * ((parti.dt + t_s)*Xi - parti.dt / (1 + d2g)) ;


    for (k=0; k<ndim;k++)
    	part[i].a_dust[k] += S * draux[k] * wkern / part[i].dt ;

    // Add Change in K.E to thermal energy generation
    if (parti.ptype == gas_type)
    	parti.dudt += S*(dvdr - 0.5*S*(1 + d2g))* wkern / parti.dt;
  }
  //-----------------------------------------------------------------------------------------------


  return;
}
*/

//=================================================================================================
//
/// Compute the value of the smoothing length of particle 'i' by iterating the relation :
/// h = h_fac*(m/rho)^(1/ndim).
/// Uses the previous value of h as a starting guess and then uses either a Newton-Rhapson solver,
/// or fixed-point iteration, to converge on the correct value of h.  The maximum tolerance used
/// for deciding whether the iteration has converged is given by the 'h_converge' parameter.
//=================================================================================================
template<int ndim, template <int> class ParticleType, class StoppingTime, class Kernel>
int DustInterpolant<ndim, ParticleType, StoppingTime, Kernel>::DoInterpolate
 (const int i,                                      ///< [in] id of particle
  const int Nneib,                                  ///< [in] No. of potential neighbours
  const FLOAT hmax,                                 ///< [in] Max. h permitted by neib list
  const std::vector<FLOAT>& m,                      ///< [in] Array of neib. masses
  const std::vector<FLOAT>& drsqd,                  ///< [in] Array of neib. distances squared
  const std::vector<typename DustInterpolant<ndim,
  	  	  	  	  	                         ParticleType,
  	  	  	  	  	                         StoppingTime,
  	  	  	  	  	                         Kernel>::DataType>& d,  ///< [in] Array of velocities
  ParticleType<ndim> &parti)                         ///< [inout] Particle i data
{
  int j;                               // Neighbour id
  int k;                               // Dimension counter
  int iteration = 0 ;
  int iteration_max = 100 ;
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT h_lower_bound = (FLOAT) 0.0;   // Lower bound on h
  FLOAT h_upper_bound = hmax;          // Upper bound on h
  FLOAT invhsqd;                       // (1 / h)^2
  FLOAT ssqd;                          // Kernel parameter squared, (r/h)^2
  FLOAT w ;                            // Kernel value
  FLOAT n ;                            // Kernel normalization
  FLOAT grho ;                         // Gas Density
  FLOAT gsound ;                       // Gas Sound Speed
  FLOAT dv[ndim] ;                     // Velocity difference
  FLOAT da[ndim] ;                     // Acceleration difference


  // Some basic sanity-checking in case of invalid input into routine
  assert(Nneib > 0);
  assert(hmax > (FLOAT) 0.0);
  assert(parti.itype != dead);

  FLOAT h = parti.h_dust ;

  // Use a guess to prevent against unknown initial smoothing length
  if (h == 0)
	  h = 0.5 * hmax ;

  // Main smoothing length iteration loop
  //===============================================================================================
  do {

    // Initialise all variables for this value of h
    iteration++;
    n = 0 ;
    grho = 0 ;
    gsound = 0 ;
    for (k=0; k<ndim; k++)
    	dv[k] = da[k] = 0 ;

    FLOAT invh = 1./h ;
    FLOAT hfactor = pow(invh, ndim) ;
    FLOAT invhsqd = invh*invh ;

    // Loop over all nearest neighbours in list to calculate density, omega and zeta.
    //---------------------------------------------------------------------------------------------
    for (j=0; j<Nneib; j++) {
      w = kern.w0_s2(drsqd[j]*invhsqd);
      n      += w ;
      grho   += m[j]*w ;
      gsound += m[j]*w*d[j].cs ;
      for (k=0; k < ndim; k++)
      {
    	  dv[k] += m[j]*w*d[j].v[k] ;
    	  da[k] += m[j]*w*d[j].a[k] ;
      }
    }
    //---------------------------------------------------------------------------------------------

   n   	*= hfactor;
   grho *= hfactor;

   FLOAT invrho = 1 / grho ;

    gsound *= invrho * hfactor ;
    for (k=0; k < ndim; k++)
     {
    	dv[k] = parti.v[k] - dv[k]*invrho*hfactor ;
    	da[k] = parti.a[k] - da[k]*invrho*hfactor ;
     }

    // If h changes below some fixed tolerance, exit iteration loop
    if (n > (FLOAT) 0.0 && h > h_lower_bound &&
        fabs(h - h_fac*pow(n,-1./ndim)) < h_converge) break;

    // Use fixed-point iteration, i.e. h_new = h_fac*(m/rho_old)^(1/ndim), for now.  If this does
    // not converge in a reasonable number of iterations (iteration_max), then assume something is
    // wrong and switch to a bisection method, which should be guaranteed to converge, albeit much
    // more slowly.  (N.B. will implement Newton-Raphson soon)
    //---------------------------------------------------------------------------------------------
    if (iteration < iteration_max) {
    	if (n > 0)
           h = h_fac*pow(n,-1./ndim);
    	else if (h < hmax)
    		h = hmax ;
    	else
    		return 0 ;
    }
    else if (iteration == iteration_max) {
      h = (FLOAT) 0.5*(h_lower_bound + h_upper_bound);
    }
    else if (iteration < 5*iteration_max) {
      if (n*pow(h,ndim) > pow(h_fac,ndim)) {
        h_upper_bound = h;
      }
      else {
        h_lower_bound = h;
      }
      h = (FLOAT) 0.5*(h_lower_bound + h_upper_bound);
    }

    if (iteration > 5*(iteration_max-2)){
      cout << "H ITERATION (DUST): " << iteration << "    h : " << h
           << "   n : " << n << "   h_upper " << h_upper_bound << "    hmax :  " << hmax
           << "   h_lower : " << h_lower_bound << "    " << hfactor << "    m : " << parti.m
           << "     " << parti.hfactor*kern.w0(0.0) << "    " << Nneib << endl;

      if (iteration == 5*iteration_max){
        string message = "Problem with convergence of h-rho iteration";
        ExceptionHandler::getIstance().raise(message);
      }
    }
    // If the smoothing length is too large for the neighbour list, exit routine and flag neighbour
    // list error in order to generate a larger neighbour list (not properly implemented yet).
    if (h > hmax) return 0;

  } while (h > h_lower_bound && h < h_upper_bound);
  //===============================================================================================

  assert(!(isinf(h)) && !(isnan(h)));
  assert(h >= h_lower_bound);


  parti.h_dust = h ;
  parti.sound = gsound ;
  parti.div_v = sqrt(DotProduct(dv,dv, ndim)) / parti.h ;

  assert(parti.h_dust > 0) ;
  assert(parti.sound > 0)  ;
  assert(parti.div_v >= 0) ;

  // Scale these so h rather than h dust can be used in the time-step calculation
  parti.div_v *= parti.h / parti.h_dust ;
  parti.sound *= parti.h / parti.h_dust ;



  //===============================================================================================
  // Compute the drag acceleration
  FLOAT t_s = t_stop(grho, 0, gsound) ;
  if (t_s == 0){
	  cout << iteration << " " << hmax << endl ;
	  cout << grho << " " << gsound << " " << h << " " << n << " " << endl ;
	  for (k=0;k<ndim;k++) cout << parti.r[k] << " " ;
	  cout <<endl;
	  for (k=0;k<ndim;k++) cout << parti.a_dust[k] << " " ;
	  cout <<endl;
	  //for (int kk(0); kk < Nneib; kk++)
	 //	  cout << "\t" << d[kk].cs << " "<< m[kk] <<" " << drsqd[kk] <<  endl ;
  }

  assert(t_s != 0) ;

  FLOAT Xi, Lambda ;
  FLOAT tau = parti.dt / t_s ;
  if (tau > 1e-3) {
    Xi      = (1 - exp(- tau)) / parti.dt ;
    Lambda  = (parti.dt + t_s) * Xi - 1;
  }
  else {
	Xi = (1 - 0.5 * tau * (1 - tau/3.)) ;
	Lambda = (1 + tau) * Xi - 1;
	Xi /= t_s ;
  }
 // Xi = 1/t_s ;
 // Lambda = 0 ;

  for (k=0; k<ndim;k++)
	parti.a_dust[k] = - dv[k] * Xi  + da[k] * Lambda ;

  // If h is invalid (i.e. larger than maximum h), then return error code (0)
  if (parti.h <= hmax) return 1;
  else return -1;
}



//=================================================================================================
//  Funtion _DustFactoryKern
/// \brief   DustFactory function definition for selecting the kernel template.
/// \details Sets up the drag force object from the parameters.
/// \author  R. A. Booth
/// \date    17/10/2015
//=================================================================================================
template<int ndim, template<int> class ParticleType, class StoppingTime, class Kernel>
class _DustFactoryKern
{
public:
DustBase<ndim>* ProcessParameters(Parameters* simparams,
							      TreeBase<ndim>* t, TreeBase<ndim>* ghost, TreeBase<ndim>* mpi_tree)
{
	map<string, int> &intparams = simparams->intparams;
	map<string, double> &floatparams = simparams->floatparams;
	map<string, string> &stringparams = simparams->stringparams;
	string KernelName = stringparams["kernel"];
	string DustForces = stringparams["dust_forces"];

	double K_D  = floatparams["drag_coeff"] ;

	if (DustForces == "test_particle")
	{
		typedef  DustTestParticle<ndim, ParticleType, StoppingTime, Kernel> dust ;
		typename dust::DI interp(StoppingTime(K_D), Kernel(KernelName),
				                 floatparams["h_fac"], floatparams["h_converge"]) ;

		DustSphNgbFinder<ndim, ParticleType>* d =
	    		new DustTestParticle<ndim,ParticleType, StoppingTime, Kernel>(interp, t, ghost) ;
#ifdef MPI_PARALLEL
		d->set_mpi_tree(mpi_tree) ;
#endif
	    return d ;
	}
	else
	{
	    string message = "Invalid option for the Dust force parameter: " +
	      stringparams["dust_forces"];
	    ExceptionHandler::getIstance().raise(message);
	}
	return NULL ;
}
};

//=================================================================================================
//  Function DustFactoryStop
/// \brief   DustFactory function definition for selecting the stopping time function
/// \details Sets up the drag force object from the parameters.
/// \author  R. A. Booth
/// \date    17/10/2015
//=================================================================================================
template<int ndim, template<int> class ParticleType, class StoppingTime>
class _DustFactoryStop
{
public:

DustBase<ndim>* ProcessParameters(Parameters * simparams,
							      TreeBase<ndim>* t, TreeBase<ndim>* ghost, TreeBase<ndim>* mpi_tree)
{
	map<string, int> &intparams = simparams->intparams;
	map<string, double> &floatparams = simparams->floatparams;
	map<string, string> &stringparams = simparams->stringparams;
	string KernelName = stringparams["kernel"];

	if (intparams["tabulated_kernel"] == 1) {
		_DustFactoryKern<ndim, ParticleType, StoppingTime, TabulatedKernel<ndim> > DF;
		return DF.ProcessParameters(simparams, t, ghost, mpi_tree);
	  }

	  else if (intparams["tabulated_kernel"] == 0) {
	    // Depending on the kernel, instantiate a different GradSph object
	    if (KernelName == "m4") {
			_DustFactoryKern<ndim, ParticleType, StoppingTime, M4Kernel<ndim> > DF;
			return DF.ProcessParameters(simparams, t, ghost, mpi_tree);
	    }
	    else if (KernelName == "quintic") {
			_DustFactoryKern<ndim, ParticleType, StoppingTime, QuinticKernel<ndim> > DF;
			return DF.ProcessParameters(simparams, t, ghost, mpi_tree);
	    }
	    else if (KernelName == "gaussian") {
			_DustFactoryKern<ndim, ParticleType, StoppingTime, GaussianKernel<ndim> > DF;
			return DF.ProcessParameters(simparams, t, ghost, mpi_tree);
	    }
	    else {
	      string message = "Unrecognised parameter : kernel = " + simparams->stringparams["kernel"];
	      ExceptionHandler::getIstance().raise(message);
	    }
	  }

	  else {
	    string message = "Invalid option for the tabulated_kernel parameter: " +
	      stringparams["tabulated_kernel"];
	    ExceptionHandler::getIstance().raise(message);
	  }
	return NULL ;
}
};


//=================================================================================================
//  Funtion DustFactory::Process Parameters
/// \brief   DustFactory function definition
/// \details Selects the appropriate
/// \author  R. A. Booth
/// \date    17/10/2015
//=================================================================================================
template<int ndim, template<int> class ParticleType>
DustBase<ndim>* DustFactory<ndim, ParticleType>::ProcessParameters
(Parameters * simparams,
TreeBase<ndim>* t,
TreeBase<ndim>* ghost,
TreeBase<ndim>* mpi_tree)
{
	map<string, int> &intparams = simparams->intparams;
	map<string, double> &floatparams = simparams->floatparams;
	map<string, string> &stringparams = simparams->stringparams;
	string DragLaw = stringparams["drag_law"];

	if (stringparams["dust_forces"] == "none")
		return NULL ;

	if (stringparams["neib_searhch"] == "bruteforce"){
	  ExceptionHandler::getIstance().raise("Dust forces are not compatible with brute force "
			                               "neighbour finding") ;
	}


	// Depending on the kernel, instantiate a different GradSph object
	if (DragLaw == "fixed") {
		_DustFactoryStop<ndim, ParticleType, FixedDrag> DF ;
		return DF.ProcessParameters(simparams, t, ghost, mpi_tree) ;
	}
	else if (DragLaw == "density") {
		_DustFactoryStop<ndim, ParticleType, DensityDrag> DF ;
		return DF.ProcessParameters(simparams, t, ghost, mpi_tree) ;
	}
	else if (DragLaw == "epstein") {
		_DustFactoryStop<ndim, ParticleType, EpsteinDrag> DF ;
		return DF.ProcessParameters(simparams, t, ghost, mpi_tree) ;	}

	else {
		string message = "Unrecognised parameter : drag_law = " + simparams->stringparams["drag_law"];
		ExceptionHandler::getIstance().raise(message);
	}
	return NULL ;
}

template class DustFactory<1, GradhSphParticle> ;
template class DustFactory<2, GradhSphParticle> ;
template class DustFactory<3, GradhSphParticle> ;
