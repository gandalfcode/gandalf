#include "Integration.h"
#include "MeshlessFV.h"

//=================================================================================================
//  MfvIntegration<ndim>::Timestep
/// Compute timestep for particle based on Courant and acceleration conditions.
//=================================================================================================
template <int ndim, template <int> class ParticleType>
DOUBLE MfvIntegration<ndim,ParticleType>::Timestep(Particle<ndim>& _part, Hydrodynamics<ndim>* hydro)
{
  MeshlessFV<ndim>* mfv = reinterpret_cast<MeshlessFV<ndim>*>(hydro) ;

  MeshlessFVParticle<ndim>& part = *reinterpret_cast<MeshlessFVParticle<ndim>*>(&_part);

  FLOAT dt_cfl = 2*courant_mult*part.h/part.vsig_max;
  FLOAT dt_grav = accel_mult*
    sqrt(part.h/sqrt(DotProduct(part.a0, part.a0, ndim) + small_number));

 if (visc_coeff > 0) {
   FLOAT dt_visc = visc_mult * 0.5*part.h*part.h/visc_coeff;
   dt_cfl = min(dt_cfl, dt_visc);
 }

  if (dt_cfl < 1e-10 || dt_grav < 1e-10) {
    cout << part.iorig << " " << part.ptype << " "
         << part.h << ", " << part.sound << " " << part.vsig_max
         << ", " << DotProduct(part.a0, part.a0, ndim) << "\n";

  }

  if (mfv->hydro_forces) return min(dt_cfl, dt_grav);
  else if (mfv->self_gravity) return dt_grav;
  else return big_number;
}


//=================================================================================================
//  MeshlessFV<ndim>::IntegrateParticles
/// Calculate or reset all quantities for all particles that reach the end of their timesteps.
//=================================================================================================
template <int ndim, template <int> class ParticleType>
void MfvIntegration<ndim,ParticleType>::AdvanceParticles
 (const int level_step,                ///< [in] Block timestep level for lowest step
  const int n,                         ///< [in] Integer time in block time struct
  const FLOAT t,                       ///< [in] Current simulation time
  const FLOAT timestep,                ///< [in] Base timestep value
  Hydrodynamics<ndim>* hydro)          ///< [inout] Pointer to Hydrodynamics object
{
  debug2("[MfvIntegration::IntegrateParticles]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("MFV_ADVANCE_PARTICLES");


  MeshlessFV<ndim>* mfv = reinterpret_cast<MeshlessFV<ndim>*>(hydro);
  MeshlessFVParticle<ndim>* partdata = mfv->GetMeshlessFVParticleArray();

  int nvar   = MeshlessFV<ndim>::nvar ;
  int irho   = MeshlessFV<ndim>::irho ;
  int ietot  = MeshlessFV<ndim>::ietot;


  // Integrate all conserved variables to end of timestep
  //-----------------------------------------------------------------------------------------------
#pragma omp parallel for default(none) shared(partdata, mfv, nvar, irho, ietot, cout)
  for (int i=0; i<mfv->Nhydro; i++) {
    MeshlessFVParticle<ndim> &part = partdata[i];
    if (part.flags.is_dead()) continue;

    // Compute time since beginning of current step
    const int nstep = pow(2, level_step - part.level);
    int dn = n%nstep;
    if (dn == 0) dn = nstep;
    const FLOAT dt = timestep*(FLOAT) dn;
    FLOAT Qcons[nvar];

    if (dn == nstep) {
      part.flags.set(active);
      for (int k=0; k<nvar; k++) Qcons[k] = part.Qcons0[k] + part.dQ[k];
    }
    else {
      part.flags.unset(active);
      for (int k=0; k<nvar; k++) Qcons[k] = part.Qcons0[k] + part.dQdt[k]*dt;
    }

    // Add the acceleration / change in energy due to gravity
    for (int k=0; k<ndim; k++) {
      Qcons[ietot] += 0.5 * dt *
          (part.a0[k]*(part.Qcons0[k] + 0.5*part.Qcons0[irho]*part.a0[k]*part.dt) +
           part.a0[k]*(Qcons[k]       + 0.5*     Qcons [irho]*part.a0[k]*part.dt));

      Qcons[k] += 0.5*(part.Qcons0[irho] + Qcons[irho])*part.a0[k]*dt;
    }

    // Add any cooling
    Qcons[ietot] -= part.cooling*dt;


    // Some sanity-checking
    //assert(isnormal(Qcons[irho]));
    //assert(isnormal(Qcons[ipress]));


    // Compute primitive values and update all main array quantities
    mfv->UpdateArrayVariables(part, Qcons);
    mfv->ComputeThermalProperties(part);
    mfv->UpdatePrimitiveVector(part);

    //---------------------------------------------------------------------------------------------
    if (!staticParticles) {
      part.flags.set(update_density);

      //-------------------------------------------------------------------------------------------
      for (int k=0; k<ndim; k++)
        part.r[k] = part.r0[k] + (FLOAT) 0.5*(part.v0[k] + part.v[k])*dt;
    }
  }

  return;
}



//=================================================================================================
//  MeshlessFV<ndim>::EndTimestep
/// Calculate or reset all quantities for all particles that reach the end of their timesteps.
//=================================================================================================
template <int ndim, template <int> class ParticleType>
void MfvIntegration<ndim,ParticleType>::EndTimestep
 (const int level_step,                ///< [in] Block timestep level for lowest step
  const int n,                         ///< [in] Integer time in block time struct
  const FLOAT t,                       ///< [in] Current simulation time
  const FLOAT timestep,                ///< [in] Base timestep value
  Hydrodynamics<ndim> *hydro)          ///< [inout] Pointer to Hydrodynamics object
{
  debug2("[MfvIntegration::EndTimestep]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("MFV_END_TIMESTEP");

  MeshlessFV<ndim>* mfv = reinterpret_cast<MeshlessFV<ndim>*>(hydro);
  MeshlessFVParticle<ndim>* partdata = mfv->GetMeshlessFVParticleArray();

  int nvar  = MeshlessFV<ndim>::nvar;
  int irho  = MeshlessFV<ndim>::irho;
  int ietot = MeshlessFV<ndim>::ietot;


  //-----------------------------------------------------------------------------------------------
#pragma omp parallel for default(none) shared(partdata, mfv, nvar, irho, ietot, cout)
  for (int i=0; i<mfv->Nhydro; i++) {
    MeshlessFVParticle<ndim> &part = partdata[i];    // Local reference to particle
    if (part.flags.is_dead()) continue;

    // Compute time since beginning of current step
    const int nstep = pow(2, level_step - part.level);
    int dn = n%nstep;
    if (dn == 0) dn = nstep;

    // If particle is at the end of its timestep
    //---------------------------------------------------------------------------------------------
    if (dn == nstep || part.flags.check(end_timestep)) {

      // Integrate all conserved quantities to end of the step (adding sums from neighbours)
      FLOAT Qcons[nvar];
      for (int var=0; var<nvar; var++) {
        Qcons[var]     = part.Qcons0[var] + part.dQ[var];
        part.dQ[var]   = (FLOAT) 0.0;
        part.dQdt[var] = (FLOAT) 0.0;
      }

      // Further update conserved quantities if computing gravitational/nbody  contributions
      for (int k=0; k<ndim; k++) {
        Qcons[ietot] += 0.5 * part.dt *
            (part.a0[k]*(part.Qcons0[k] + 0.5*part.Qcons0[irho]*part.a0[k]*part.dt) +
             part.a [k]*(Qcons[k]       + 0.5*     Qcons [irho]*part.a [k]*part.dt));

        Qcons[k] += 0.5 * part.dt * (part.Qcons0[irho]*part.a0[k] + Qcons[irho]*part.a[k]);
      }
      Qcons[ietot] += 0.5 *
         (DotProduct(part.a0, part.rdmdt, ndim) + DotProduct(part.a, part.rdmdt, ndim));

      // Add any cooling:
      Qcons[ietot] -= part.cooling * part.dt;

      // Compute primitive values and update all main array quantities
      mfv->UpdateArrayVariables(part, Qcons);
      mfv->ComputeThermalProperties(part);
      mfv->UpdatePrimitiveVector(part);

      // Update all values to the beginning of the next step
      for (int k=0; k<ndim; k++) part.r0[k]     = part.r[k];
      for (int k=0; k<ndim; k++) part.v0[k]     = part.v[k];
      for (int k=0; k<ndim; k++) part.a0[k]     = part.a[k];
      for (int k=0; k<ndim; k++) part.rdmdt0[k] = part.rdmdt[k];
      for (int k=0; k<ndim; k++) part.rdmdt[k]  = 0.0;
      for (int k=0; k<nvar; k++) part.Qcons0[k] = Qcons[k];
      part.cooling = 0;
      part.dt      = part.dt_next;
      part.dt_next = 0;
      part.flags.set(active);
      part.flags.unset(end_timestep);

    }
    //---------------------------------------------------------------------------------------------
    else {
      part.flags.unset(active);
    }
    //---------------------------------------------------------------------------------------------

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  MeshlessFV<ndim>::EndTimestep
/// Calculate or reset all quantities for all particles that reach the end of their timesteps.
//=================================================================================================
template <int ndim, template <int> class ParticleType>
int MfvIntegration<ndim,ParticleType>::CheckTimesteps
 (const int level_diff_max,            ///< [in] Max. allowed SPH neib dt diff
  const int level_step,                ///< [in] Level of base timestep
  const int n,                         ///< [in] Integer time in block time struct
  const FLOAT timestep,                ///< [in] Timestep
  Hydrodynamics<ndim>* hydro)          ///< [inout] Pointer to Hydrodynamics object
{
  int activecount = 0;                 // No. of newly active particles
  int nvar = MeshlessFV<ndim>::nvar;
  double tstep = timestep;

  MeshlessFV<ndim>* mfv = reinterpret_cast<MeshlessFV<ndim>*>(hydro);
  MeshlessFVParticle<ndim> *mfvdata = mfv->GetMeshlessFVParticleArray();

  debug2("[MeshlessFV::CheckTimesteps]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("MESHLESS_CHECK_TIMESTEPS");

  if (mfv->timestep_limiter != "simple") return 0;

  //-----------------------------------------------------------------------------------------------
  #pragma omp parallel for default(none) shared(cout, mfv, mfvdata, nvar, tstep) reduction(+:activecount)
  for (int i=0; i<mfv->Nhydro; i++) {
    MeshlessFVParticle<ndim>& part = mfvdata[i];
    if (part.flags.is_dead()) continue;

    // Compute time since beginning of current step
    const int nstep = pow(2, level_step - part.level);
    int dn = n%nstep;
    if (dn == 0) dn = nstep;

    // Check if neighbour timesteps are too small.  If so, then reduce timestep if possible
    if (part.levelneib - part.level > level_diff_max) {
      const int level_new = part.levelneib - level_diff_max;
      const int nnewstep  = pow(2, level_step - level_new);

      // Force recalculation of fluxes for particles at the end of their step
      if (dn%nnewstep == 0 && dn != nstep) {
        part.level = level_new;
        part.dt = (FLOAT) nstep * tstep;

        // Use current predicted value for dQ
        for (int var=0; var<nvar; var++) {
          part.dQ[var] = part.dt * part.dQdt[var];
        }

        part.flags.set(active);
        part.flags.set(sm_limiter);
      }
    }
  }
  //-----------------------------------------------------------------------------------------------

  return activecount;
}



//=================================================================================================
//  SphLeapfrogKDK::SetActiveParticles
/// Set or unset the active flag for all particles based upon whther the particles need a force
/// calculation this timestep.
//=================================================================================================
template <int ndim, template <int> class ParticleType>
void MfvIntegration<ndim, ParticleType>::SetActiveParticles
 (const int level_step,                ///< [in] Block timestep level for lowest step
  const int n,                         ///< [in] Current timestep number
  Hydrodynamics<ndim>* hydro)          ///< [inout] Pointer to Hydrodynamics object
{
  MeshlessFV<ndim>* mfv = reinterpret_cast<MeshlessFV<ndim>*>(hydro);
  MeshlessFVParticle<ndim> *mfvdata = mfv->GetMeshlessFVParticleArray() ;

#pragma omp parallel for default(none) shared(mfvdata, mfv)
  for (int i=0; i<mfv->Nhydro; i++) {
    Particle<ndim>& part = mfvdata[i];

    // Compute time since beginning of current step
    const int nstep = pow(2, level_step - part.level);
    int dn = n%nstep;
    if (dn == 0) dn = nstep;

    // Force calculation is at end of step
    if (part.flags.check(end_timestep) || dn == nstep || dn == 0) {
      part.flags.set(active);
    }
    else {
      part.flags.unset(active);
    }
  }
}

template class MfvIntegration<1,MeshlessFVParticle>;
template class MfvIntegration<2,MeshlessFVParticle>;
template class MfvIntegration<3,MeshlessFVParticle>;
