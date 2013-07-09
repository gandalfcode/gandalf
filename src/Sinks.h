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

  StarParticle<ndim> *star;         ///< Pointer to connected star particle
  int istar;                        ///< i.d. of connected star particle
  int Ngas;                         ///< No. of gas particles inside sink
  DOUBLE radius;                    ///< Softening/sink radius of particle
  DOUBLE dt;                        ///< Particle timestep
  DOUBLE dmdt;                      ///< Accretion rate
  DOUBLE menc;                      ///< ..
  DOUBLE mmax;                      ///< ..
  DOUBLE racc;                      ///< Accretion radius
  DOUBLE ketot;                     ///< Internal kinetic energy
  DOUBLE gpetot;                    ///< Internal grav. pot. energy
  DOUBLE rotketot;                  ///< Internal rotational kinetic energy
  DOUBLE utot;                      ///< Internal energy accrted by sink
  DOUBLE taccrete;                  ///< ..
  DOUBLE trad;                      ///< ..
  DOUBLE trot;                      ///< ..
  DOUBLE tvisc;                     ///< ..
  DOUBLE angmom[3];                 ///< Internal sink angular momentum


  // Star particle constructor to initialise all values
  // --------------------------------------------------------------------------
  SinkParticle()
  {
    radius = 0.0;
    dt = 0.0;
    dmdt = 0.0;
    racc = 0.0;
    ketot = 0.0;
    gpetot = 0.0;
    rotketot = 0.0;
    utot = 0.0;
    for (int k=0; k<3; k++) angmom[k] = 0.0;
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

  // Function prototypes
  // --------------------------------------------------------------------------
  void AllocateMemory(int);
  void DeallocateMemory(void);
  void SearchForNewSinkParticles(int, Sph<ndim> *, Nbody<ndim> *);
  void CreateNewSinkParticle(int, Sph<ndim> *, Nbody<ndim> *);
  void AccreteMassToSinks(Sph<ndim> *, Nbody<ndim> *, int, DOUBLE, int);
  //void UpdateSystemProperties(void);

  // Local class variables
  // --------------------------------------------------------------------------
  bool allocated_memory;            ///< Has sink memory been allocated?
  int Nsink;                        ///< No. of sink particles
  int Nsinkmax;                     ///< Max. no. of sink particles
  int sink_particles;               ///< Using sink particles?
  int create_sinks;                 ///< Create new sink particles?
  int smooth_accretion;             ///< Use smooth accretion?
  FLOAT alpha_ss;                   ///< Shakura-Sunyaev alpha viscosity
  FLOAT rho_sink;                   ///< Sink formation density
  FLOAT sink_radius;                ///< New sink radius (in units of h)
  FLOAT smooth_accrete_frac;        ///< ..
  FLOAT smooth_accrete_dt;          ///< ..
  string sink_radius_mode;          ///< Sink radius mode

  SinkParticle<ndim> *sink;         ///< Main sink particle array

};
#endif
