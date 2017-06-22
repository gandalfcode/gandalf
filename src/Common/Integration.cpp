#include "Integration.h"
#include "Debug.h"


//=================================================================================================
//  TimeIntegration::CheckBoundaries
/// Check all particles to see if any have crossed the simulation bounding box.
/// If so, then move the particles to their new location on the other side of the periodic box.
//=================================================================================================
template <int ndim>
void TimeIntegration<ndim>::CheckBoundaries
 (DomainBox<ndim>& simbox,             ///< Domain box object
  Hydrodynamics<ndim>* hydro)          ///< Pointer to SPH object
{
  // If all boundaries are open, immediately return to main loop
  if (simbox.boundary_lhs[0] == openBoundary && simbox.boundary_rhs[0] == openBoundary &&
      simbox.boundary_lhs[1] == openBoundary && simbox.boundary_rhs[1] == openBoundary &&
      simbox.boundary_lhs[2] == openBoundary && simbox.boundary_rhs[2] == openBoundary) return;

  debug2("[TimeIntegration::CheckBoundaries]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("CHECK_BOUNDARIES");

  // Loop over all particles and check if any lie outside the periodic box.
  // If so, then re-position with periodic wrapping.
  //===============================================================================================
#pragma omp parallel for default(none) shared(simbox,hydro)
  for (int i=0; i<hydro->Nhydro; i++) {
    Particle<ndim>& part = hydro->GetParticlePointer(i);

    // --------------------------------------------------------------------------------------------
    for (int k=0; k<ndim; k++) {

      // Check if particle has crossed LHS boundary
      //-------------------------------------------------------------------------------------------
      if (part.r[k] < simbox.min[k]) {

        // Check if periodic boundary
        if (simbox.boundary_lhs[k] == periodicBoundary) {
          part.r[k]  += simbox.size[k];
          part.r0[k] += simbox.size[k];
        }

        // Check if wall or mirror boundary
        if (simbox.boundary_lhs[k] == mirrorBoundary || simbox.boundary_lhs[k] == wallBoundary) {
          part.r[k]  = (FLOAT) 2.0*simbox.min[k] - part.r[k];
          part.r0[k] = (FLOAT) 2.0*simbox.min[k] - part.r0[k];
          part.v[k]  = -part.v[k];
          part.v0[k] = -part.v0[k];
          part.a[k]  = -part.a[k];
          part.a0[k] = -part.a0[k];
        }

      }

      // Check if particle has crossed RHS boundary
      //-------------------------------------------------------------------------------------------
      if (part.r[k] > simbox.max[k]) {

        // Check if periodic boundary
        if (simbox.boundary_rhs[k] == periodicBoundary) {
          part.r[k]  -= simbox.size[k];
          part.r0[k] -= simbox.size[k];
        }

        // Check if wall or mirror boundary
        if (simbox.boundary_rhs[k] == mirrorBoundary || simbox.boundary_rhs[k] == wallBoundary) {
          part.r[k]  = (FLOAT) 2.0*simbox.max[k] - part.r[k];
          part.r0[k] = (FLOAT) 2.0*simbox.max[k] - part.r0[k];
          part.v[k]  = -part.v[k];
          part.v0[k] = -part.v0[k];
          part.a[k]  = -part.a[k];
          part.a0[k] = -part.a0[k];
        }

      }


    }
    //---------------------------------------------------------------------------------------------

  }
  //===============================================================================================

  return;
}



template class TimeIntegration<1>;
template class TimeIntegration<2>;
template class TimeIntegration<3>;
