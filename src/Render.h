// ============================================================================
// Render.h
// ============================================================================


#ifndef _RENDER_H_
#define _RENDER_H_


#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <math.h>
#include "SphParticle.h"
#include "Sph.h"
#include "SphSnapshot.h"
#include "Exception.h"
#include "Parameters.h"
#include "SphKernel.h"
#include "InlineFuncs.h"
#include "Debug.h"
using namespace std;


class RenderBase
{
public:
  static RenderBase* RenderFactory(int ndim, SimulationBase* sim);

  virtual int CreateColumnRenderingGrid(int, int, string, string, string, string,
                float, float, float, float, float* values,
                int Ngrid, SphSnapshotBase &, float& scaling_factor)=0;
  virtual int CreateSliceRenderingGrid(int, int, string, string, string, string, string,
      float, float, float, float, float, float* values,
      int Ngrid,
      SphSnapshotBase &, float& scaling_factor)=0;
};


// ============================================================================
// Class Render
// ============================================================================
template <int ndim>
class Render : public RenderBase
{
 public:

  // Constructor and Destructor
  // --------------------------------------------------------------------------
  Render(SimulationBase* sim);
  ~Render();

  // Subroutine prototypes
  // --------------------------------------------------------------------------
  int CreateColumnRenderingGrid(int, int, string, string, string, string,
				float, float, float, float, float* values, 
				int Ngrid, SphSnapshotBase &, float& scaling_factor);
  int CreateSliceRenderingGrid(int, int, string, string, string, string, string,
			       float, float, float, float, float, float* values, 
			       int Ngrid,
			       SphSnapshotBase &, float& scaling_factor);


  // ..
  // --------------------------------------------------------------------------
  Sph<ndim>* sph;


};


#endif
