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
#include "NeighbourManager.h"

//=================================================================================================
//  Class DustSphNgbFinder
/// \brief   DustSphNgbFinder class definition.
/// \details Template base class that does neighbour searching for drag forces. Two methods are
///          provided: (1) FindNeigbAndDoInterp,
///                    (2) FindNeigbAndDoForces,
///          (1) finds the gas neighbours to dust particles for interpolating gas properties to
///          the location of dust particles, while (2) finds the dust/gas neighbours for pair-wise
///          force calculations in the full two-fluid limit that can include the back reaction on
///          the gas.
///
///          The implementer should provide a template interpolation object if (1) is used, which
///          must have define the nested type DataType and method DoInterpolate. DataType must be
///          constructible from ParticleType and should collect the interpolation data needed.
///          DoInterpolation should compute the dust smoothing length h_dust, and compute the drag
///          forces from the interpolated data. See DustInterpolant for an example.
///
///          If (2) is used, then the implementer must provide the ComputeDragForces function that
///          computes the drag forces for a given particle using the neighbours found by
///          FindNeibAndDoForces. See DustSemiImplictForces for an example.
/// \author  R. A. Booth
/// \date    17/10/2015
//=================================================================================================
template<int ndim, template<int> class ParticleType>
class DustSphNgbFinder
: public DustBase<ndim>
{
	using DustBase<ndim>::timing ;
	vector<NeighbourManager<ndim,ParticleType<ndim> > > neibmanagerbuf;
public:
	DustSphNgbFinder(TreeBase<ndim> * t,TreeBase<ndim> * gt=NULL)
	: _tree(t), _ghosttree(gt)
	{
	} ;
	virtual ~DustSphNgbFinder() { } ;

#ifdef MPI_PARALLEL
	void set_mpi_tree(TreeBase<ndim>* t)
	{ mpighosttree = t ; }
#endif

protected:
	template<class Interp>
	void FindNeibAndDoInterp(int ,int, ParticleType<ndim> *,Typemask, Interp&) ;

	template<class ForceCalc>
	void FindNeibAndDoForces(int ,int, ParticleType<ndim> *,
				             const ParticleTypeRegister&, ForceCalc&) ;
private:

	TreeBase<ndim>* _tree, *_ghosttree ;   ///< Pointer to neighbour tree
#if defined MPI_PARALLEL
	TreeBase<ndim>* mpighosttree;        ///< Pointer to pruned tree arrays
#endif
protected:

};


//=================================================================================================
//  Class DustInterpolation
/// \brief   DustInterpolation class definition.
/// \details Interpolater class used for computing drag forces in the test particle limit. Provides
///          a data type that can be used for collecting the necessary density, sound speed
///          velocity and acceleration data from a gas particle, along with the interpolation and
///          semi-implicit force calculation in the test particle limit.
/// \author  R. A. Booth
/// \date    17/10/2015
//=================================================================================================
template<int ndim, template <int> class ParticleType, class StoppingTime, class Kernel>
class DustInterpolant
{
  //===============================================================================================
  //  Class DustTestPartInterp
  /// \brief   DustTestPartInterp class definition.
  /// \details Collector class that grabs the required interpolation data for test particle drag
  ///          force calculation.
  /// \author  R. A. Booth
  /// \date    17/10/2015
  //===============================================================================================
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
//  Class DustSemiImplictForces
/// \brief   DustSemiImplictForces class definition.
/// \details Class that computes the drag forces between gas and dust particles in the full two-
///          fluid approximation including the back-reaction of the dust on the gas from the given
///          neighbour list.
/// \author  R. A. Booth
/// \date    5/11/2015
//=================================================================================================
template<int ndim, template <int> class ParticleType, class StoppingTime, class Kernel>
class DustSemiImplictForces
{
public:

  DustSemiImplictForces(const StoppingTime& ts,const Kernel &k, bool UseEnergyTerm)
   : kern(k), t_stop(ts), _use_energy_term(UseEnergyTerm)
  { } ;

  void ComputeDragForces(const int, const int, const std::vector<int>&,
		                     ParticleType<ndim>&, const std::vector<ParticleType<ndim> >&,
		                     std::vector<FLOAT>&) ;


private:
  Kernel kern ;
  StoppingTime t_stop ;
  bool _use_energy_term ;
};


//=================================================================================================
//  Class DustFull
/// \brief   DustFull class definition.
/// \details Computes the drag force on the dust and it's back-reaction on the gas using the
///          semi-implicit update.
/// \author  R. A. Booth
/// \date    17/10/2015
//=================================================================================================
template<int ndim, template<int> class ParticleType, class StoppingTime, class Kernel>
class DustFull
: public DustSphNgbFinder<ndim,ParticleType>
{
  using DustSphNgbFinder<ndim,ParticleType>::FindNeibAndDoForces ;
public:
  typedef DustSemiImplictForces<ndim, ParticleType, StoppingTime, Kernel>  DF ;

  DustFull(DF Forces, const ParticleTypeRegister& types,
		  TreeBase<ndim> * t, TreeBase<ndim> * gt=NULL)
    :  DustSphNgbFinder<ndim, ParticleType>(t, gt),
       _types(types),
	   _Forces(Forces)
  { } ;

  void UpdateAllDragForces(int NPart, int Ntot, Particle<ndim> *sph_gen){

    debug2("[DustFull::UpdateAllDragForces]") ;

    ParticleType<ndim>* sphdata = reinterpret_cast<ParticleType<ndim>*>(sph_gen) ;

    FindNeibAndDoForces(NPart, Ntot, sphdata, _types, _Forces) ;
  }
private:
  ParticleTypeRegister _types ;
  DF _Forces ;
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

		Typemask mask ;
		mask[gas_type] = true ;

		FindNeibAndDoInterp(NPart, Ntot, sphdata, mask, _interp) ;
	}
private:
  DI _interp ;
};


//=================================================================================================
//  Member function DustSphNgbFinder::FindNgbInterpolate
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
  vector<TreeCellBase<ndim> > celllist;    // List of active cells

#ifdef MPI_PARALLEL
  double twork = timing->RunningTime();  // Start time (for load balancing)
#endif

  debug2("[DustSphNgbFinder::FindNeibAndDoInterp]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("DUST_GAS_INTERPOLATE_FORCES");

  // Find list of all cells that contain active particles
  cactive = _tree->ComputeActiveCellList(celllist);
  assert(cactive <= _tree->MaxNumCells());


  // Set-up all OMP threads
  //===============================================================================================
#pragma omp parallel default(none) shared(cactive,celllist,cout,sphdata, mask, Interp)
 {
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
      TreeCellBase<ndim>& cell = celllist[cc];
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
          if (sphdata[j].flags.is_dead()){
        	Nneib-- ;
        	continue ;
          }

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
          if (activepart[j].ptype != dust_type) continue ;

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
#ifdef MPI_PARALLEL
 twork = timing->RunningTime() - twork;
 int Nactivetot=0;
 _tree->AddWorkCost(celllist, twork, Nactivetot) ;
#ifdef OUTPUT_ALL
 cout << "Time computing dust pair-wise forces: " << twork << "     Nactivetot : " << Nactivetot << endl;
#endif
#endif

  return;
}


//=================================================================================================
//  Member function DustSphNgbFinder::FindNgbForces
// \brief    Driver routine for pairwise particle forces
/// \details This function finds the pairs of particles that are involved in a given force
///          calculation.
/// \author  R. A. Booth
/// \date    20/10/2015
//=================================================================================================
template<int ndim, template<int> class ParticleType>
template<class ForceCalc>
void DustSphNgbFinder<ndim, ParticleType>::FindNeibAndDoForces
(int Nhydro,                              ///< [in] No. of SPH particles
 int Ntot,                                ///< [in] No. of SPH + ghost particles
 ParticleType<ndim> *sphdata,             ///< [inout] Pointer to SPH ptcl array
 const ParticleTypeRegister& types,       ///< [in] Type data for particles
 ForceCalc& Forces)                       ///< [in] Force calculation functor
{
  using std::vector ;

  int cactive;                             // No. of active cells
  vector<TreeCellBase<ndim> > celllist;    // List of active cells


#ifdef MPI_PARALLEL
  double twork = timing->RunningTime();  // Start time (for load balancing)
#endif

#ifdef _OPENMP
  int Nthreads  = omp_get_max_threads() ;
#else
  int Nthreads  = 1 ;
#endif
  for (int t = neibmanagerbuf.size(); t < Nthreads; ++t)
    neibmanagerbuf.push_back(NeighbourManager<ndim,
                                              ParticleType<ndim> >(types,_tree->MaxKernelRange(),
                                                                   _tree->GetDomain()));


  debug2("[DustSphNgbFinder::FindNeibForces]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("DUST_GAS_PAIRWISE_FORCES");


  // Find list of all cells that contain active particles
  cactive = _tree->ComputeActiveCellList(celllist);

  // If there are no active cells, return to main loop
  if (cactive == 0) {
    return;
  }

  _tree->UpdateAllHmaxValues(sphdata);

  // Set-up all OMP threads
  //===============================================================================================
  vector<FLOAT> a_drag(ndim*Ntot) ;            // temporary to hold the drag accelerations

#pragma omp parallel default(none) shared(cactive,celllist,sphdata,types,Forces, Nhydro, Ntot, a_drag)
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
    vector<int>                 activelist(_tree->MaxNumPartInLeafCell()); // Ids of Active parts
    vector<ParticleType<ndim> > activepart(_tree->MaxNumPartInLeafCell()); // Local array of parts
    vector<int>                 levelneib(Ntot,0);                         // Ngb t-step level
    NeighbourManager<ndim,ParticleType<ndim> >& neibmanager = neibmanagerbuf[ithread];


    // Loop over all active cells
    //=============================================================================================
#pragma omp for schedule(guided)
    for (cc=0; cc<cactive; cc++) {
      TreeCellBase<ndim>& cell = celllist[cc];

      // Find list of active particles in current cell
      Nactive = _tree->ComputeActiveParticleList(cell,sphdata, &(activelist[0]));

      // Make local copies of active particles
      for (j=0; j<Nactive; j++) {
        activepart[j] = sphdata[activelist[j]];
      }


      // Compute neighbour list for cell from real and periodic ghost particles
      neibmanager.clear();
      _tree->ComputeNeighbourAndGhostList(cell, neibmanager);
#ifdef MPI_PARALLEL
        // Ghosts are already in the mpi tree
        mpighosttree->ComputeNeighbourList(cell, neibmanager);
#endif
      neibmanager.EndSearch(cell,sphdata);

      const int Nneib_cell = neibmanager.GetNumAllNeib();

      // Loop over all active particles in the cell
      //-------------------------------------------------------------------------------------------
      for (j=0; j<Nactive; j++) {
        i = activelist[j];

        if (types[activepart[j].ptype].drag_forces) {

          Typemask dragmask ;
          for (k=0; k<Ntypes; k++)  dragmask[k] = types[activepart[j].ptype].dragmask[k];

          const bool do_pair_once=false;
          int* sphlist_temp;
          ParticleType<ndim>* neibpart_temp;

          const int Nneib=neibmanager.GetParticleNeib(activepart[j],dragmask,&sphlist_temp,&neibpart_temp,do_pair_once);

          vector<int> sphlist(sphlist_temp,sphlist_temp+Nneib);
          vector<ParticleType<ndim> > neibpart (neibpart_temp,neibpart_temp+Nneib_cell);
      	  // Compute all neighbour contributions to hydro forces
      	  vector<FLOAT> acc(ndim) ;
          Forces.ComputeDragForces(i,Nneib,sphlist,activepart[j],neibpart, acc);
          for (k=0; k<ndim; k++)
            a_drag[i*ndim + k] = acc[k] ;

        }
      }

      // Add all active particles contributions to main array
      for (j=0; j<Nactive; j++) {
    	i = activelist[j];
    	sphdata[i].dudt = activepart[j].dudt ;
    	sphdata[i].div_v = activepart[j].div_v ;
    	sphdata[i].sound = activepart[j].sound ;
      }
    }
    //===============================================================================================

#pragma omp for schedule(guided)
    for (i=0; i<Nhydro; i++) {
      sphdata[i].levelneib = max(sphdata[i].levelneib, levelneib[i]);
      for (k=0; k<ndim; k++)
        sphdata[i].a[k] += a_drag[i*ndim + k] ;
    }
  }

  // Compute time spent in routine and in each cell for load balancing
#ifdef MPI_PARALLEL
  twork = timing->RunningTime() - twork;
  int Nactivetot=0;
  _tree->AddWorkCost(celllist, twork, Nactivetot) ;
#ifdef OUTPUT_ALL
  cout << "Time computing dust pair-wise forces: " << twork << "     Nactivetot : " << Nactivetot << endl;
#endif
#endif

  return;
}


//=================================================================================================
//   Member function DustInterpolant:DoInterpolate
/// \brief   Interpolation routine for calculating dust smoothing length and drag forces.
/// \details Compute the smoothing length from the number density of gas neighbours using
///          fixed-point iteration, in a similar way to the SPH density calculation. Once the
///          local gas properties have been determined the stopping-time is calculated and the drag
///          acceleration is averaged over the time-step. The gas sound speed and dust-gas relative
///          velocity are then stored in sound and div_v for the time-step calculation.
/// \author  R. A. Booth
/// \date    20/10/2015
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
  FLOAT h_lower_bound = (FLOAT) 0.0;   // Lower bound on h
  FLOAT h_upper_bound = hmax;          // Upper bound on h
  FLOAT w ;                            // Kernel value
  FLOAT n ;                            // Kernel normalization
  FLOAT grho ;                         // Gas Density
  FLOAT gsound ;                       // Gas Sound Speed
  FLOAT dv[ndim] ;                     // Velocity difference
  FLOAT da[ndim] ;                     // Acceleration difference


  // Some basic sanity-checking in case of invalid input into routine
  assert(Nneib > 0);
  assert(hmax > (FLOAT) 0.0);
  assert(!parti.flags.is_dead());

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

  // Predict the relative velocity
  for (k=0; k<ndim; k++)
    dv[k] += da[k] * parti.dt ;

  //===============================================================================================
  // Compute the drag acceleration
  FLOAT t_s = t_stop(grho, 0, gsound) ;

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

  for (k=0; k<ndim; k++)
    parti.a[k] += - dv[k] * Xi  + da[k] * Lambda ;

  // If h is invalid (i.e. larger than maximum h), then return error code (0)
  if (parti.h <= hmax) return 1;
  else return -1;
}


//=================================================================================================
//   Member function DustSemiImplictForces::ComputeDragForces
/// \brief   Pairwise force calculation for semi-implicit drag forces.
/// \details Compute the drag acceleration on a particles due to its neighbours using a pair-wise
///          force projected along the line joining the particles in an angular-momentum conserving
///          way. The drag acceleration is averaged over the particle's time-step. If the particle
///          is a dust particle the maximum sound-speed of its neighbours and dust-gas relative
///          velocity are stored for the time-step computation. For gas particles, the change in
///          kinetic energy is added to the internal energy to conserve energy.
/// \author  R. A. Booth
/// \date    6/11/2015
//=================================================================================================
template<int ndim, template <int> class ParticleType, class StoppingTime, class Kernel>
void DustSemiImplictForces<ndim, ParticleType, StoppingTime, Kernel>::ComputeDragForces
(const int i,                                       ///< [in] id of particle
 const int Nneib,                                   ///< [in] No. of neibs in neibpart array
 const std::vector<int>& neiblist,                  ///< [in] id of gather neibs in neibpart
 ParticleType<ndim>& parti,                         ///< [inout] Particle i data
 const std::vector<ParticleType<ndim> >& neibpart,  ///< [in] Neighbour particle data
 std::vector<FLOAT>& a_drag)                        ///< [out] drag acceleration
{
  int j;                               // Neighbour list id
  int jj;                              // Aux. neighbour counter
  int k;                               // Dimension counter
  FLOAT draux[ndim];                   // Relative position vector
  FLOAT dv[ndim];                      // Relative velocity vector
  FLOAT da[ndim];                      // Relative acceleration vector
  FLOAT dvdr;                          // Dot product of dv and dr
  FLOAT dadr;                          // Dot product of da and dr
  FLOAT wkern;                         // Value of drag kernel function for part i
  FLOAT S ;                            // Drag term

  // Some basic sanity-checking in case of invalid input into routine
  assert(!parti.flags.is_dead());

  if (parti.ptype == dust_type){
	  parti.sound = 0;
	  parti.div_v = 0;
  }

  FLOAT invh_i =  1/parti.h ;

  // Loop over all potential neighbours in the list
  //-----------------------------------------------------------------------------------------------
  for (jj=0; jj<Nneib; jj++) {
    j = neiblist[jj];
    assert(!neibpart[j].flags.is_dead());

    FLOAT invh_j =  1/neibpart[j].h ;

    for (k=0; k<ndim; k++) draux[k] = neibpart[j].r[k] - parti.r[k];
    const FLOAT drmag =  sqrt(DotProduct(draux,draux,ndim));


    if (parti.ptype == gas_type)
      wkern = pow(invh_i, ndim)*kern.wdrag(drmag*invh_i);
    else
      wkern = pow(invh_j, ndim)*kern.wdrag(drmag*invh_j);

    wkern *= neibpart[j].m / neibpart[j].rho ;

    for (k=0; k<ndim; k++) {
    	if (drmag>0) draux[k] /= drmag;
    	dv[k] = neibpart[j].v[k] - parti.v[k] ;
    	da[k] = neibpart[j].a[k] - parti.a[k] ;
    }
    dvdr = DotProduct(draux, dv, ndim);
    dadr = DotProduct(draux, da, ndim);

    //===============================================================================================
    // Compute stopping time
    FLOAT gsound, grho, drho;
    if (parti.ptype == gas_type){
      gsound = parti.sound ;
      grho = parti.rho ;
      drho = neibpart[j].rho ;
    } else {
      gsound = neibpart[j].sound ;
      grho = neibpart[j].rho ;
      drho = parti.rho ;
      parti.sound = max(parti.sound, gsound) ;
      parti.div_v = max(parti.div_v, sqrt(DotProduct(dv,dv,ndim))/parti.h);
    }


    FLOAT t_s = t_stop(grho, drho, gsound) ;
    assert(t_s > 0) ;

    //---------------------------------------------------------------------------------------------
    // Evaluate the drag term
    FLOAT rho = drho + grho ;
    FLOAT tau = parti.dt / t_s ;
    FLOAT Xi, Lambda ;
    if (tau > 1e-3) {
      Xi = (1 - exp(- tau)) / (parti.dt * rho) ;
      Lambda = (parti.dt + t_s)*Xi - 1 / rho ;
    } else {
      Xi = (1 - 0.5 * tau * (1 - tau/3.)) / rho ;
      Lambda = (1 + tau) * Xi - 1 / rho;
      Xi /= t_s ;
    }

    // Predict the relative velocity
    dvdr += parti.dt * dadr ;

    S = (dvdr * Xi - dadr * Lambda) ;

    for (k=0; k<ndim; k++)
    	a_drag[k] += ndim * neibpart[j].rho * S * draux[k] * wkern ;

    // Add Change in K.E to thermal energy generation
    if (_use_energy_term && parti.ptype == gas_type)
    	parti.dudt += ndim * neibpart[j].rho * S *(dvdr - 0.5*rho*S*parti.dt)* wkern ;
  }
  //-----------------------------------------------------------------------------------------------

  return;
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
DustBase<ndim>* ProcessParameters(Parameters* simparams, ParticleTypeRegister& types,
							      TreeBase<ndim>* t, TreeBase<ndim>* ghost, TreeBase<ndim>* mpi_tree)
{
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

		DustSphNgbFinder<ndim, ParticleType>* d = new dust(interp, t, ghost) ;
#ifdef MPI_PARALLEL
		d->set_mpi_tree(mpi_tree) ;
#endif
	    return d ;
	}
	else if (DustForces == "full_twofluid") {
	  typedef DustFull<ndim, ParticleType, StoppingTime, Kernel> dust ;

	  StoppingTime t_s(K_D) ; Kernel kern(KernelName) ;
	  bool IntegrateEnergy = stringparams["eos"] != "isothermal" ;

	  typename dust::DF Forces(t_s, kern, IntegrateEnergy) ;
	  DustSphNgbFinder<ndim, ParticleType>* d =  new dust(Forces, types, t, ghost) ;
#ifdef MPI_PARALLEL
	  d->set_mpi_tree(mpi_tree) ;
#endif
	  return d ;
	}
	else {
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

DustBase<ndim>* ProcessParameters(Parameters * simparams, ParticleTypeRegister& types,
							      TreeBase<ndim>* t, TreeBase<ndim>* ghost, TreeBase<ndim>* mpi_tree)
{
	map<string, int> &intparams = simparams->intparams;
	map<string, string> &stringparams = simparams->stringparams;
	string KernelName = stringparams["kernel"];

	if (intparams["tabulated_kernel"] == 1) {
		_DustFactoryKern<ndim, ParticleType, StoppingTime, TabulatedKernel<ndim> > DF;
		return DF.ProcessParameters(simparams, types, t, ghost, mpi_tree);
	  }

	  else if (intparams["tabulated_kernel"] == 0) {
	    // Depending on the kernel, instantiate a different GradSph object
	    if (KernelName == "m4") {
			_DustFactoryKern<ndim, ParticleType, StoppingTime, M4Kernel<ndim> > DF;
			return DF.ProcessParameters(simparams, types, t, ghost, mpi_tree);
	    }
	    else if (KernelName == "quintic") {
			_DustFactoryKern<ndim, ParticleType, StoppingTime, QuinticKernel<ndim> > DF;
			return DF.ProcessParameters(simparams, types, t, ghost, mpi_tree);
	    }
	    else if (KernelName == "gaussian") {
			_DustFactoryKern<ndim, ParticleType, StoppingTime, GaussianKernel<ndim> > DF;
			return DF.ProcessParameters(simparams, types, t, ghost, mpi_tree);
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
//  Function DustFactory::Process Parameters
/// \brief   DustFactory function definition
/// \details Selects the appropriate
/// \author  R. A. Booth
/// \date    17/10/2015
//=================================================================================================
template<int ndim, template<int> class ParticleType>
DustBase<ndim>* DustFactory<ndim, ParticleType>::ProcessParameters
(Parameters * simparams,
CodeTiming * timing,
ParticleTypeRegister& types,
TreeBase<ndim>* t,
TreeBase<ndim>* ghost,
TreeBase<ndim>* mpi_tree)
{
	map<string, int> &intparams = simparams->intparams;
	map<string, string> &stringparams = simparams->stringparams;
	string DragLaw = stringparams["drag_law"];

	if (stringparams["dust_forces"] == "none")
		return NULL ;

	if (stringparams["neib_searhch"] == "bruteforce"){
	  ExceptionHandler::getIstance().raise("Dust forces are not compatible with brute force "
			                               "neighbour finding") ;
	}

	if (intparams["dimensionless"] == 0){
	  ExceptionHandler::getIstance().raise("Error: Non-dimensionless simulations with dust are "
			  	  	  	  	  	  	  	   "not currently supported") ;
	}

	DustBase<ndim> * dust_forces ;

	// Depending on the kernel, instantiate a different GradSph object
	if (DragLaw == "fixed") {
		_DustFactoryStop<ndim, ParticleType, FixedDrag> DF ;
		dust_forces = DF.ProcessParameters(simparams, types, t, ghost, mpi_tree) ;
	}
	else if (DragLaw == "density") {
		_DustFactoryStop<ndim, ParticleType, DensityDrag> DF ;
		dust_forces = DF.ProcessParameters(simparams, types, t, ghost, mpi_tree) ;
	}
	else if (DragLaw == "epstein") {
		_DustFactoryStop<ndim, ParticleType, EpsteinDrag> DF ;
		dust_forces = DF.ProcessParameters(simparams, types, t, ghost, mpi_tree) ;
	}
	else if (DragLaw == "LP2012") {
	    _DustFactoryStop<ndim, ParticleType, LP12_Drag> DF ;
	    dust_forces = DF.ProcessParameters(simparams, types, t, ghost, mpi_tree) ;
	}
	else {
		string message = "Unrecognised parameter : drag_law = " + simparams->stringparams["drag_law"];
		ExceptionHandler::getIstance().raise(message);
		return NULL ;
	}

	dust_forces->timing = timing ;
	return dust_forces ;
}

template class DustFactory<1, GradhSphParticle> ;
template class DustFactory<2, GradhSphParticle> ;
template class DustFactory<3, GradhSphParticle> ;
