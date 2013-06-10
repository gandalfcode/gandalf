//=============================================================================
//  Sinks.h
//  Main sink particle class
//=============================================================================


#ifndef _SINKS_H_
#define _SINKS_H_


#include <string>
#include "Precision.h"
#include "Constants.h"
#include "Parameters.h"
#include "SphKernel.h"
#include "NbodyParticle.h"
#include "StarParticle.h"
#include "SystemParticle.h"
#include "SphParticle.h"
#include "Sph.h"
#include "Nbody.h"
using namespace std;



//=============================================================================
//  Class SinkParticle
/// \brief   Individual sink particle data structure
/// \details Main parent N-body particle data structure.  All main other 
///          N-body particle types (e.g. stars, systems, sinks) are derived 
///          from this class.
/// \author  D. A. Hubber
/// \date    15/04/2013
//=============================================================================
template <int ndim>
class SinkParticle
{
 public:

  bool active;                      ///< Flag if active (i.e. recompute step)
  int level;                        ///< Current timestep level    
  int nstep;                        ///< Integer step-size of particle
  DOUBLE r[ndim];                   ///< Position
  DOUBLE v[ndim];                   ///< Velocity
  DOUBLE a[ndim];                   ///< Acceleration
  DOUBLE adot[ndim];                ///< Time derivative of acceleration (jerk)
  DOUBLE m;                         ///< Star mass
  DOUBLE h;                         ///< Smoothing length
  DOUBLE invh;                      ///< 1 / h
  DOUBLE radius;                    ///< Softening/sink radius of particle
  DOUBLE dt;                        ///< Particle timestep
  DOUBLE dmdt;                      ///< Accretion rate
  DOUBLE racc;                      ///< Accretion radius
  DOUBLE ketot;                     ///< Internal kinetic energy
  DOUBLE gpetot;                    ///< Internal grav. pot. energy
  DOUBLE rotketot;                  ///< Internal rotational kinetic energy


  // Star particle constructor to initialise all values
  // --------------------------------------------------------------------------
  SinkParticle()
  {
    active = false;
    level = 0;
    nstep = 0;
    for (int k=0; k<ndim; k++) r[k] = 0.0;
    for (int k=0; k<ndim; k++) v[k] = 0.0;
    for (int k=0; k<ndim; k++) a[k] = 0.0;
    for (int k=0; k<ndim; k++) adot[k] = 0.0;
    m = 0;
    h = 0;
    invh = 0.0;
    radius = 0.0;
    dt = 0.0;
    dmdt = 0.0;
    racc = 0.0;
    ketot = 0.0;
    gpetot = 0.0;
    rotketot = 0.0;
  } 

};



//=============================================================================
//  Class Sinks
/// \brief   Main sink particle class.
/// \details Main sink particle class for searching for and creating new 
///          sinks, and for controlling the accretion of SPH particles 
///          onto sink particles.
/// \author  D. A. Hubber
/// \date    08/06/2013
//=============================================================================
template<int ndim>
class Sinks
{
 public:

  Sinks();
  ~Sinks();

  void AllocateMemory(void);
  void DeallocateMemory(void);
  void SearchForNewSinkParticles(int, Sph<ndim> *, Nbody<ndim> *);
  void CreateNewSinkParticle(int);
  void AccreteMassToSinks(void);
  void UpdateSinkProperties(void);
  void UpdateStarProperties(void);
  //void UpdateSystemProperties(void);

  
  int Nsink;                        ///< No. of sink particles
  int Nsinkmax;                     ///< Max. no. of sink particles
  FLOAT alpha_ss;                   ///< Shakura-Sunyaev alpha viscosity
  FLOAT rho_sink;                   ///< Sink formation density
  FLOAT sink_radius;                ///< New sink radius (in units of h)

  SinkParticle<ndim> *sink;         ///< Main sink particle array

};
#endif
