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
#include <limits>
#include <vector>
#include "CodeTiming.h"
#include "Debug.h"
#include "DragLaws.h"
#include "Dust.h"
#include "MeshlessFV.h"
#ifdef MPI_PARALLEL
#include "MpiControl.h"
#endif
#include "Particle.h"
#include "Precision.h"
#include "Tree.h"
#include "NeighbourManager.h"

/* Velocity at the start of the kick */
template<int ndim>
double get_initial_velocity(const GradhSphParticle<ndim>& part, int k) {
  return part.v[k] - 0.5*part.dt * part.a0[k];
}

template<int ndim>
double get_initial_velocity(const MeshlessFVParticle<ndim>& part, int k) {
  return part.v0[k] ;
}

/* Difference in velocity at the start of the kick */
template<int ndim, class NeibType>
double get_velocity_difference(const GradhSphParticle<ndim>& part,
    const NeibType& neib, int k) {
  return (part.v[k] - neib.v[k]) - 0.5*part.dt *(part.a0[k] - neib.a0[k]);
}

template<int ndim, class NeibType>
double get_velocity_difference(const MeshlessFVParticle<ndim>& part,
    const NeibType& neib, int k) {
  return part.v0[k] - neib.v0[k] ;
}

/* Acceleration during the kick */
template<int ndim>
double get_total_accel(const GradhSphParticle<ndim>& part, int k) {
  return part.a[k] ;
}

/* For the meshless, use the time-averaged momentum change */
template<int ndim>
double get_total_accel(const MeshlessFVParticle<ndim>& part, int k) {
  double a_grav  = 0.5*(part.a[k] + part.a0[k]*part.Qcons0[MeshlessFV<ndim>::irho]/part.m);
  double a_hydro = 0 ;
  if (part.dt > 0)
    a_hydro = part.dQ[k] / (part.m * part.dt) ;
  return a_grav + a_hydro ;
}

/* Acceleration at the start of the time-step (might be used to correct the initial velocity */
template<int ndim>
double get_old_accel(const GradhSphParticle<ndim>& part, int k) {
  return part.a0[k] ;
}

/* Not used for the meshless */
template<int ndim>
double get_old_accel(const MeshlessFVParticle<ndim>& part, int k) {
  return 0 ;
}

template<int ndim>
void update_particle(GradhSphParticle<ndim>& part, FLOAT* acc, FLOAT dudt, double) {
  for(int k=0; k<ndim; k++)
    part.a[k] += acc[k] ;

  if(part.ptype == gas_type)
    part.dudt += dudt ;
}

template<int ndim>
void update_particle(MeshlessFVParticle<ndim>& part, FLOAT* acc, FLOAT dudt, double dt) {

  // dudt gets updated using particles exported to other processors
  dudt += part.dudt ;

  if (part.ptype == gas_type) {
    // Add change in Kinetic energy to dudt to get total change in total energy

    for (int k=0; k<ndim; k++) {
      double v0 = get_initial_velocity(part,k) + get_total_accel(part, k)*dt;
      dudt += acc[k] * (v0 + acc[k]*dt/2) ;
    }

    part.dQ[MeshlessFV<ndim>::ietot]   += part.m * dudt * dt;
    part.dQdt[MeshlessFV<ndim>::ietot] += part.m * dudt ;
  }
  else {
    assert(dudt == 0) ;
    part.vsig_max = part.h*std::abs(part.div_v) + part.sound ;
  }

  for(int k=0; k<ndim; k++) {
    part.dQ[k]   += part.m * acc[k] * dt ;
    part.dQdt[k] += part.m * acc[k] ;
  }

  // Update the particle position according to the drag terms.
  for(int k=0; k<ndim; k++) {
    part.v[k] = part.v0[k] + part.dQ[k] + part.a0[k]*dt ;
    part.r[k] = part.r0[k] + 0.5*(part.v0[k] + part.v[k])*dt;
    part.flags.set(update_density) ;
  }

  // Throw away now we have added it to the conserved arrays
  part.dudt = 0 ;
}



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
#ifdef MPI_PARALLEL
  using DustBase<ndim>::mpicontrol ;
#endif
  vector<NeighbourManager<ndim,ParticleType<ndim> > > neibmanagerbuf;
public:
  DustSphNgbFinder(Parameters* params, TreeBase<ndim> * t,TreeBase<ndim> * gt=NULL)
: _tree(t), _ghosttree(gt), w_tstep(std::numeric_limits<DOUBLE>::signaling_NaN())
{
    std::string simtype = params->stringparams["sim"] ;
    if (simtype.find("sph") != std::string::npos) {
      if (params->stringparams["sph_integration"] == "lfkdk")
        w_tstep = 0.5;
      if (params->stringparams["sph_integration"] == "lfdkd")
        w_tstep = 1.0;
    }
    else if (simtype == "meshlessfv" || simtype == "mfvmuscl") {
      w_tstep = 1.0 ;
    }

    // Check that we set w_tstep
    if (w_tstep != w_tstep) {
      std::string message =
          "Problem with in DustSphNgbFinder constructor, time-integration not recognised";
      ExceptionHandler::getIstance().raise(message);
    }
} ;
  virtual ~DustSphNgbFinder() { } ;

#ifdef MPI_PARALLEL
  void set_mpi_tree(TreeBase<ndim>* t)
  { mpighosttree = t ; }
#endif

protected:
  template<class Interp>
  void FindNeibAndDoInterp(Hydrodynamics<ndim>*, Typemask, Interp&) ;

  template<class ForceCalc>
  void FindNeibAndDoForces(Hydrodynamics<ndim>*,
      const ParticleTypeRegister&, ForceCalc&) ;
private:
  TreeBase<ndim>* _tree, *_ghosttree ;   ///< Pointer to neighbour tree
#if defined MPI_PARALLEL
  TreeBase<ndim>* mpighosttree;          ///< Pointer to pruned tree arrays
#endif

  DOUBLE w_tstep;                        ///< Weight for current / previous time-step
protected:
  DOUBLE drag_timestep(ParticleType<ndim>& part) {
    return w_tstep*part.dt + (1 - w_tstep)*part.dt_next ;
  }

private:
  std::vector<std::vector<FLOAT> > dudt_buf ;
};


//===============================================================================================
//  Class DustTestPartInterp
/// \brief   DustTestPartInterp class definition.
/// \details Collector class that grabs the required interpolation data for test particle drag
///          force calculation.
/// \author  R. A. Booth
/// \date    17/10/2015
//===============================================================================================
template<int ndim, template <int> class ParticleType>
struct DustTestPartInterp
{
  DustTestPartInterp() : cs(0), m(0), ptype(0), hrangesqd(0) { } ;
  DustTestPartInterp(const ParticleType<ndim>& p){
    cs = p.sound ;
    m  = p.m ;
    ptype = p.ptype;
    flags = p.flags;
    hrangesqd = p.hrangesqd;
    for (int k=0; k < ndim; ++k){
      r[k] = p.r[k] ;
      v[k]  = p.v[k] ;
      v0[k] = p.v0[k];
      a[k]  = get_total_accel(p, k) ;
      a0[k] = get_old_accel(p, k) ;
    }
  }

  FLOAT cs, m, hrangesqd;
  FLOAT r[ndim];
  FLOAT v[ndim] ;
  FLOAT v0[ndim];
  FLOAT a[ndim] ;
  FLOAT a0[ndim];
  int ptype ;
  type_flag flags;

  typedef ParticleType<ndim> BaseParticle;

  static const int NDIM = ndim ;
};

template<int ndim>
void reflect(DustTestPartInterp<ndim,GradhSphParticle>&, int, double)
{
  ExceptionHandler::getIstance().raise("You should not use mirror boundaries with dust!!!!");
}

template<int ndim>
void reflect(DustTestPartInterp<ndim,MeshlessFVParticle>&, int, double)
{
  ExceptionHandler::getIstance().raise("You should not use mirror boundaries with dust!!!!");
}


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

public:

  DustInterpolant(const StoppingTime& ts,const Kernel &k, FLOAT h_factor, FLOAT h_conv )
: kern(k), t_stop(ts), h_fac(h_factor), h_converge(h_conv)
{ } ;

  typedef DustTestPartInterp<ndim,ParticleType> DataType;

  int DoInterpolate(ParticleType<ndim>&, const NeighbourList<DataType>&, DOUBLE, FLOAT) ;

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

  void ComputeDragForces(ParticleType<ndim>&, NeighbourList<ParticleType<ndim> >&, DOUBLE,
      FLOAT *) ;

  bool NeedEnergyUpdate() const { return _use_energy_term ; }

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

  DustFull(Parameters* params,DF Forces, const ParticleTypeRegister& types,
      TreeBase<ndim> * t, TreeBase<ndim> * gt=NULL)
  :  DustSphNgbFinder<ndim, ParticleType>(params,t, gt),
     _types(types),
     _Forces(Forces)
     { } ;

  DF _Forces ;
  void UpdateAllDragForces(Hydrodynamics<ndim>* hydro){

    debug2("[DustFull::UpdateAllDragForces]") ;

    FindNeibAndDoForces(hydro, _types, _Forces) ;
  }
private:
  ParticleTypeRegister _types ;
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

  DustTestParticle(Parameters* params, DI Interp, TreeBase<ndim> * t, TreeBase<ndim> * gt=NULL)
  :  DustSphNgbFinder<ndim, ParticleType>(params, t, gt),
     _interp(Interp)
     { } ;
  void UpdateAllDragForces(Hydrodynamics<ndim>* hydro)
  {
    debug2("[DustTestParticle::UpdateAllDragForces]") ;

    Typemask mask ;
    mask[gas_type] = true ;

    FindNeibAndDoInterp(hydro, mask, _interp) ;
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
(Hydrodynamics<ndim>* hydro,                ///< [in] Hydro class
 Typemask mask,                            ///< [in] Types to include in the interpolation
 Interpolant& Interp)                      ///< [in] Interpolating function
{
  using std::vector ;

  typedef typename Interpolant::DataType InterpData;

  static vector<NeighbourManager<ndim, InterpData> > neibbuf;

  debug2("[DustSphNgbFinder::FindNeibAndDoInterp]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("DUST_GAS_INTERPOLATE_FORCES");

  ParticleType<ndim>* sphdata = hydro->template GetParticleArray<ParticleType>();

  int cactive;                             // No. of active tree cells
  vector<TreeCellBase<ndim> > celllist;    // List of active cells

#ifdef MPI_PARALLEL
  double twork = timing->RunningTime();  // Start time (for load balancing)
#endif

#ifdef _OPENMP
  int Nthreads  = omp_get_max_threads() ;
#else
  int Nthreads  = 1 ;
#endif
  for (int t = neibbuf.size(); t < Nthreads; ++t)
    neibbuf.push_back(NeighbourManager<ndim,InterpData>(hydro->types,_tree->MaxKernelRange(),
        _tree->GetDomain()));

  // Find list of all cells that contain active particles
  cactive = _tree->ComputeActiveCellList(celllist);
  assert(cactive <= _tree->MaxNumCells());


  // Set-up all OMP threads
  //===============================================================================================
#pragma omp parallel default(none) shared(cactive,celllist,cout,sphdata,mask,Interp,hydro,neibbuf)
  {
#if defined _OPENMP
      const int ithread = omp_get_thread_num();
#else
      const int ithread = 0;
#endif
    int celldone;                              // Flag if cell is done
    int cc;                                    // Aux. cell counter
    int j;                                     // Aux. particle counter
    int Nactive;                               // No. of active particles in cell
    int okflag;                                // Flag if particle is done
    FLOAT hmax;                                // Maximum smoothing length
    FLOAT hrangesqd;                           // Smoothing range  (squared)
    vector<int> activelist(_tree->MaxNumPartInLeafCell()) ;  // Local array of active particles ids
    vector<ParticleType<ndim> > activepart(_tree->MaxNumPartInLeafCell()) ;  // Local array of active particles
    NeighbourManager<ndim,InterpData>& neibmanager = neibbuf[ithread];

    FLOAT kernrangesqd = pow(_tree->MaxKernelRange(),2) ;


    // Loop over all active cells
    //=============================================================================================
#pragma omp for schedule(guided)
    for (cc=0; cc<cactive; cc++) {
      TreeCellBase<ndim> cell = celllist[cc];
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
        cell.hmax = hmax;


        // Compute neighbour list for cell from particles on all trees
        neibmanager.set_target_cell(cell);
        _tree->ComputeGatherNeighbourList(cell,sphdata,hmax,neibmanager);
        _ghosttree->ComputeGatherNeighbourList(cell,sphdata,hmax,neibmanager);
#ifdef MPI_PARALLEL
        mpighosttree->ComputeGatherNeighbourList(cell,sphdata,hmax,neibmanager);
#endif
        neibmanager.EndSearchGather(cell, sphdata);


        // Loop over all active particles in the cell
        //-----------------------------------------------------------------------------------------
        for (j=0; j<Nactive; j++) {

          // Skip non-dust particles
          if (activepart[j].ptype != dust_type) continue ;

          // Set gather range as current h multiplied by some tolerance factor
          hrangesqd = kernrangesqd*hmax*hmax;

          NeighbourList<InterpData> neiblist =
              neibmanager.GetParticleNeibGather(activepart[j],mask,hrangesqd);

#if defined(VERIFY_ALL)
          neibmanager.VerifyNeighbourList(activelist[j], hydro->Nhydro, sphdata, "gather");
          neibmanager.VerifyReducedNeighbourList(activelist[j], neiblist, hydro->Nhydro,
                                                 sphdata, mask, "gather");
#endif

//---------------------------------------------------------------------------------------

          DOUBLE dt_drag = drag_timestep(activepart[j]);

          // Compute smoothing length and other gather properties for ptcl i
          okflag = Interp.DoInterpolate(activepart[j], neiblist, dt_drag, hmax) ;

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
(Hydrodynamics<ndim>* hydro,              ///< [in] Hydro class
 const ParticleTypeRegister& types,       ///< [in] Type data for particles
 ForceCalc& Forces)                       ///< [in] Force calculation functor
{
  using std::vector ;

  ParticleType<ndim>* sphdata = hydro->template GetParticleArray<ParticleType>();

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
  for (int t = neibmanagerbuf.size(); t < Nthreads; ++t) {
    neibmanagerbuf.push_back(NeighbourManager<ndim,
        ParticleType<ndim> >(types,_tree->MaxKernelRange(),
            _tree->GetDomain()));

    dudt_buf.push_back(std::vector<FLOAT>()) ;
  }

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
  vector<FLOAT> a_drag(ndim*hydro->Ntot, 0) ;          // temporary to hold the drag accelerations
  vector<FLOAT>& dudt = dudt_buf[0];                   // Use a reference for the total heating


#ifdef MPI_PARALLEL
  // Zero the update for the thermal energy of ghosts so we can tell which ones have had updates
  // from other processors
  if (Forces.NeedEnergyUpdate())
    for (int i=hydro->Nhydro; i < hydro->Nhydro+hydro->NPeriodicGhost; i++)
      sphdata[i].dudt = 0 ;

  mpighosttree->UpdateAllHmaxValues(sphdata);
#endif

#pragma omp parallel default(none) shared(cactive,celllist,sphdata,types,Forces,hydro,a_drag, dudt, Nthreads)
  {
#if defined _OPENMP
      const int ithread = omp_get_thread_num();
#else
      const int ithread = 0;
#endif
      int cc;                                      // Aux. cell counter
      int i;                                       // Particle id
      int j;                                       // Aux. particle counter
      int Nactive;                                 // ..
      vector<int>                 activelist(_tree->MaxNumPartInLeafCell()); // Ids of Active parts
      vector<ParticleType<ndim> > activepart(_tree->MaxNumPartInLeafCell()); // Local array of parts
      NeighbourManager<ndim,ParticleType<ndim> >& neibmanager = neibmanagerbuf[ithread];
      vector<FLOAT>& dudt_local = dudt_buf[ithread];

      dudt_local.resize(hydro->Ntot) ;
      for (int i=0; i < hydro->Ntot; i++)
        dudt_local[i] = 0 ;

      // Loop over all active cells
      //=============================================================================================
#pragma omp for schedule(guided)
      for (cc=0; cc<cactive; cc++) {
        TreeCellBase<ndim>& cell = celllist[cc];

        // Find list of active particles in current cell
        Nactive = _tree->ComputeActiveParticleList(cell,sphdata, &(activelist[0]));

        // Make local copies of active particles
        for (j=0; j<Nactive; j++)
          activepart[j] = sphdata[activelist[j]];

        // Compute neighbour list for cell from real and periodic ghost particles
        neibmanager.set_target_cell(cell);
        _tree->ComputeNeighbourAndGhostList(cell, neibmanager);
#ifdef MPI_PARALLEL
        // Ghosts are already in the mpi tree
        mpighosttree->ComputeNeighbourList(cell, neibmanager);
#endif
        neibmanager.EndSearch(cell,sphdata);

        // Initialize the change in energy
        for (j=0; j<Nactive; j++)
          activepart[j].dudt = 0;

        int Nneib =  neibmanager.GetNumAllNeib() ;
        for (j=0; j<Nneib; j++)
          neibmanager[j].dudt = 0;

        // Loop over all active particles in the cell
        //-------------------------------------------------------------------------------------------
        for (j=0; j<Nactive; j++) {
          i = activelist[j];

          if (types[activepart[j].ptype].drag_forces) {

            Typemask dragmask = types[activepart[j].ptype].dragmask;

            const bool do_pair_once=false;
            NeighbourList<ParticleType<ndim> > neiblist =
                neibmanager.GetParticleNeib(activepart[j],dragmask,do_pair_once);

#if defined(VERIFY_ALL)
            neibmanager.VerifyNeighbourList(i, hydro->Nhydro, sphdata, "all");
            neibmanager.VerifyReducedNeighbourList(i, neiblist, hydro->Nhydro, sphdata, dragmask, "all");
#endif

            DOUBLE dt_drag = drag_timestep(activepart[j]);

            Forces.ComputeDragForces(activepart[j],neiblist, dt_drag, &(a_drag[i*ndim]));
          }

        }

        // Add all active particles contributions to main array
        for (j=0; j<Nactive; j++) {
          i = activelist[j];
          sphdata[i].div_v = activepart[j].div_v ;
          sphdata[i].sound = activepart[j].sound ;
          dudt_local[i] += activepart[j].dudt ;
        }
        for (j=0; j<Nneib; j++) {
          std::pair<int,ParticleType<ndim>*> neighbour=neibmanager.GetNeibI(j);
          dudt_local[neighbour.first] += neighbour.second->dudt ;
        }
      }

      if (Forces.NeedEnergyUpdate()) {
        // Collect together the total thermal update
        // We skip the zeroth index as we dudt_buf[0] as the array for the sum.
#pragma omp for
        for (i=0; i<hydro->Ntot; i++) {
          for (int j=1; j < Nthreads; j++)
            dudt[i] += dudt_buf[j][i] ;
        }
      }
      //===============================================================================================

#ifdef MPI_PARALLEL
      // Communicate back the dudt contributions from particles on external processors

      if (Forces.NeedEnergyUpdate()) {
#pragma omp master
        {
          std::list<int> copy_back;
          for(int n=0; n<hydro->Nmpighost; n++) {
            int i = hydro->Nhydro + hydro->NPeriodicGhost + n ;

            if (dudt[i] > 0) {
              sphdata[i].dudt = dudt[i] ;
              copy_back.push_back(i) ;
            }
          }

          mpicontrol->UpdateMpiGhostParents(copy_back, hydro, update_dust_parents);
        }

        // Barrier here to ensure MPI is finished.
#pragma omp barrier
      }
#endif



#pragma omp for
      for (i=0; i<hydro->Nhydro; i++) {
        update_particle(sphdata[i], &a_drag[i*ndim], dudt[i], drag_timestep(sphdata[i]));
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
(ParticleType<ndim> &parti,         ///< [inout] Particle i data
 const NeighbourList<typename DustInterpolant<ndim,
                     ParticleType,
                     StoppingTime,
                     Kernel>::DataType>& ngbs,  ///< [in] Neighbours
DOUBLE dt,                         ///< [in] Time to average drag force over
FLOAT hmax)                        ///< [in] Maximum smoothing length
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
  FLOAT dr[ndim] ;                     // Position difference
  FLOAT dv[ndim] ;                     // Velocity difference
  FLOAT da[ndim] ;                     // Acceleration difference
  FLOAT ssqd;

  typedef typename DustInterpolant<ndim, ParticleType,StoppingTime,Kernel>::DataType DataType;

  // Some basic sanity-checking in case of invalid input into routine
  assert(hmax > (FLOAT) 0.0);
  assert(!parti.flags.is_dead());

  FLOAT h = parti.h_dust ;

  int Nneib = ngbs.size();

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
      const DataType &ngb = ngbs[j];
      for (k=0; k<ndim; k++) dr[k] = ngb.r[k] - parti.r[k];
      ssqd = invhsqd*DotProduct(dr, dr, ndim);
      w = kern.w0_s2(ssqd);
      n      += w ;
      grho   += ngb.m*w ;
      gsound += ngb.m*w*ngb.cs ;
      for (k=0; k < ndim; k++) {
        // Get the velocity at the start of the kick and the acceleration
        dv[k] += ngb.m*w*get_velocity_difference(parti, ngb, k) ;
        da[k] += ngb.m*w*(get_total_accel(parti, k) - ngb.a[k]) ;
      }
    }
    //---------------------------------------------------------------------------------------------

    n   	*= hfactor;
    grho *= hfactor;

    FLOAT invrho = 1 / grho ;

    gsound *= invrho * hfactor ;
    for (k=0; k < ndim; k++) {
      dv[k] *= invrho*hfactor ;
      da[k] *= invrho*hfactor ;
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

  } while (h > h_lower_bound && h <= h_upper_bound);
  //===============================================================================================

  assert(!(isinf(h)) && !(isnan(h)));
  assert(h >= h_lower_bound);


  parti.h_dust = h ;
  parti.sound = gsound ;
  parti.div_v = sqrt(DotProduct(dv,dv, ndim)) / parti.h ;

  assert(parti.h_dust > 0) ;
  assert(parti.sound > 0)  ;
  assert(parti.div_v >= 0) ;

  // Scale the sound speed so h rather than h dust can be used in the time-step calculation
  parti.div_v *= parti.h / parti.h_dust ;
  parti.sound *= parti.h / parti.h_dust ;

  // Predict the relative velocity
  for (k=0; k<ndim; k++)
    dv[k] += da[k] * dt ;

  //===============================================================================================
  // Compute the drag acceleration
  FLOAT t_s = t_stop(grho, 0, gsound) ;

  assert(t_s != 0) ;

  FLOAT Xi, Lambda ;
  FLOAT tau = dt / t_s ;
  if (tau > 1e-3) {
    Xi      = (1 - exp(- tau)) / dt ;
    Lambda  = (dt + t_s) * Xi - 1;
  }
  else {
    Xi = (1 - 0.5 * tau * (1 - tau/3.)) ;
    Lambda = (1 + tau) * Xi - 1;
    Xi /= t_s ;
  }

  FLOAT a[ndim] ;
  for (k=0; k<ndim; k++)
    a[k] = - dv[k] * Xi  + da[k] * Lambda ;

  update_particle(parti, a, 0, dt) ;

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
(ParticleType<ndim>& parti,                         ///< [inout] Particle i data
 NeighbourList<ParticleType<ndim> >& neiblist,      ///< [in] List of neighbours
 DOUBLE dt,                                          ///< [in] Time to average drag force over
 FLOAT* a_drag)                                     ///< [out] drag acceleration
 {
  int j;                               // Neighbour list id
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

  if (parti.ptype != gas_type){
    parti.sound = 0;
    parti.div_v = 0;
  }

  FLOAT invh_i =  1/parti.h ;

  double _adrag[ndim] ;
  double norm = 0;
  for(k=0; k < ndim; k++) _adrag[k] = 0;

  // Loop over all potential neighbours in the list
  //-----------------------------------------------------------------------------------------------
  int Nneib = neiblist.size();
  for (j=0; j < Nneib; j++) {
    assert(!neiblist[j].flags.is_dead());

    FLOAT invh_j =  1/neiblist[j].h ;

    for (k=0; k<ndim; k++) draux[k] = parti.r[k] - neiblist[j].r[k];
    const FLOAT drmag =  sqrt(DotProduct(draux,draux,ndim));


    if (parti.ptype != dust_type)
      wkern = pow(invh_i, ndim)*kern.wdrag(drmag*invh_i);
    else
      wkern = pow(invh_j, ndim)*kern.wdrag(drmag*invh_j);

    wkern *= neiblist[j].m / neiblist[j].rho ;
    norm  += wkern ;

    for (k=0; k<ndim; k++) {
      if (drmag>0) draux[k] /= drmag;
      dv[k] = get_velocity_difference(parti, neiblist[j], k) ;
      da[k] = get_total_accel(parti, k) - get_total_accel(neiblist[j], k) ;
    }
    dvdr = DotProduct(draux, dv, ndim);
    dadr = DotProduct(draux, da, ndim);

    //===============================================================================================
    // Compute stopping time
    FLOAT gsound, grho, drho;
    if (parti.ptype == gas_type){
      gsound = parti.sound ;
      grho = parti.rho ;
      drho = neiblist[j].rho ;
    } else {
      gsound = neiblist[j].sound ;
      grho = neiblist[j].rho ;
      drho = parti.rho ;
      parti.sound = max((FLOAT) parti.sound, gsound) ;
      parti.div_v = max(parti.div_v, sqrt(DotProduct(dv,dv,ndim))/parti.h);
    }


    FLOAT t_s = t_stop(grho, drho, gsound) ;
    assert(t_s > 0) ;

    //---------------------------------------------------------------------------------------------
    // Evaluate the drag term
    FLOAT rho = drho + grho ;
    FLOAT tau = dt / t_s ;
    FLOAT Xi, Lambda ;
    if (tau > 1e-3) {
      Xi = (1 - exp(- tau)) / (dt * rho) ;
      Lambda = (dt + t_s)*Xi - 1 / rho ;
    } else {
      Xi = (1 - 0.5 * tau * (1 - tau/3.)) / rho ;
      Lambda = (1 + tau) * Xi - 1 / rho;
      Xi /= t_s ;
    }

    // Predict the relative velocity
    dvdr += dt * dadr ;

    S = (dvdr * Xi - dadr * Lambda) ;

    for (k=0; k<ndim; k++)
      _adrag[k] -= ndim * neiblist[j].rho * S * draux[k] * wkern ;
  }
  //-----------------------------------------------------------------------------------------------

  // Save the drag acceleration
  for(k=0; k < ndim; k++) a_drag[k] = _adrag[k];

  if (_use_energy_term) {
    // Change in (specific) Kinetic energy due to drag forces
    double dEk_dt = 0 ;
    for (k=0; k<ndim; k++) {
      double v0 = get_initial_velocity(parti,k) + get_total_accel(parti, k)*dt;
      dEk_dt += a_drag[k] * (v0 + a_drag[k]*dt/2) ;
    }

    if (parti.ptype == dust_type) {
      // Spread the change in energy for the dust amongst its gas neighbours.
      for (j=0; j < Nneib; j++) {
        FLOAT invh_j =  1/neiblist[j].h ;

        const FLOAT drmag = Distance(parti.r, neiblist[j].r, ndim);

        FLOAT wkern = pow(invh_j, ndim)*kern.wdrag(drmag*invh_j);
        wkern /= norm * neiblist[j].rho ;

        neiblist[j].dudt -= parti.m * wkern * dEk_dt ;
      }
    } else {
      // Keep the change in energy for the gas particle
      parti.dudt -= dEk_dt;
    }
  }
  return;
    }




//=================================================================================================
//  Function _DustFactoryKern
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
    map<string, int> &intparams = simparams->intparams;
    string KernelName = stringparams["kernel"];
    string DustForces = stringparams["dust_forces"];

    double K_D  = floatparams["drag_coeff"] ;

    if (DustForces == "test_particle")
    {
      typedef  DustTestParticle<ndim, ParticleType, StoppingTime, Kernel> dust ;
      typename dust::DI interp(StoppingTime(K_D), Kernel(KernelName),
          floatparams["h_fac"], floatparams["h_converge"]) ;

      DustSphNgbFinder<ndim, ParticleType>* d = new dust(simparams, interp, t, ghost) ;
#ifdef MPI_PARALLEL
      d->set_mpi_tree(mpi_tree) ;
#endif
      return d ;
    }
    else if (DustForces == "full_twofluid") {

      if (intparams["Nlevels"] > 1) {
        string message = "Error: Full Two-fluid dust does not work with block timesteps." ;
        ExceptionHandler::getIstance().raise(message);
      }

      typedef DustFull<ndim, ParticleType, StoppingTime, Kernel> dust ;

      StoppingTime t_s(K_D) ; Kernel kern(KernelName) ;
      bool IntegrateEnergy = stringparams["eos"] != "isothermal" ;

      typename dust::DF Forces(t_s, kern, IntegrateEnergy) ;
      DustSphNgbFinder<ndim, ParticleType>* d =  new dust(simparams, Forces, types, t, ghost) ;
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
    DomainBox<ndim>& simbox,
    TreeBase<ndim>* t,
    TreeBase<ndim>* ghost,
    TreeBase<ndim>* mpi_tree)
    {
  map<string, int> &intparams = simparams->intparams;
  map<string, string> &stringparams = simparams->stringparams;
  string DragLaw = stringparams["drag_law"];

  if (stringparams["dust_forces"] == "none")
    return NULL ;

  if (IsAnyBoundaryReflecting(simbox)) {
    ExceptionHandler::getIstance().raise(
        "Error: Dust does not work with reflecting boundaries");
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

template class DustFactory<1, MeshlessFVParticle> ;
template class DustFactory<2, MeshlessFVParticle> ;
template class DustFactory<3, MeshlessFVParticle> ;
